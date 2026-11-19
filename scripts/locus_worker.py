#!/usr/bin/env python3
"""Per-locus worker for Nextflow FINEMAP_LOCUS with coordinated node-local /tmp staging.

Coordinates transient multi-GB LD matrix (`matrix.npy`) staging across concurrent
workers on the same compute node via an `flock`-guarded reservation table in
`$SLURM_TMPDIR` / `/tmp`, sweeps orphaned `gwasfm_*` temp directories from
cancelled jobs, spills excess concurrent workers cleanly to `scratch_base` (and
retries transparently if node-local `/tmp` fills mid-write), runs `prep_locus.py`
followed by each requested `scripts/methods/run_<method>.R` script, and always
removes the transient `matrix.npy` directory in `finally:`.
"""

from __future__ import annotations

import argparse
import errno
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_LOCUS_RESERVE_BYTES = int(5.5 * (1024 ** 3))   # 5.5 GiB for standard <=9.5 Mb locus
MHC_LOCUS_RESERVE_BYTES = int(9.0 * (1024 ** 3))       # 9.0 GiB for extended 8.5 Mb MHC locus
DEFAULT_MIN_HEADROOM_BYTES = int(4.0 * (1024 ** 3))    # Leave >=4 GiB untouched on node /tmp
LEGACY_ORPHAN_AGE_SECONDS = 120.0

ENOSPC_MARKERS = (
    "not enough free space",
    "no space left on device",
    f"[errno {errno.ENOSPC}]",
    "errno 28",
)


@dataclass
class WorkdirClaim:
    work_dir: Path
    using_node_tmp: bool
    cleanup_on_exit: bool
    lock_path: Path | None = None
    state_path: Path | None = None
    claim_id: str | None = None


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _dir_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    try:
        for root, _dirs, files in os.walk(path):
            for fname in files:
                try:
                    total += os.path.getsize(os.path.join(root, fname))
                except OSError:
                    pass
    except OSError:
        pass
    return total


def _sweep_stale_tmp_dirs(tmp_root: Path, uid: int, now: float | None = None) -> int:
    """Delete orphaned gwasfm_* directories owned by `uid` whose creator PID exited."""
    now_ts = time.time() if now is None else now
    removed = 0
    try:
        entries = list(tmp_root.glob("gwasfm_*"))
    except OSError:
        return 0

    for entry in entries:
        try:
            st = entry.lstat()
        except OSError:
            continue
        if st.st_uid != uid or not entry.is_dir() or entry.is_symlink():
            continue
        pid_file = entry / ".owner_pid"
        stale = False
        if pid_file.is_file():
            try:
                owner_pid = int(pid_file.read_text(encoding="utf-8").strip())
                stale = not _pid_is_alive(owner_pid)
            except (OSError, ValueError):
                stale = (now_ts - st.st_mtime) >= LEGACY_ORPHAN_AGE_SECONDS
        else:
            # Legacy directory created before .owner_pid tracking
            stale = (now_ts - st.st_mtime) >= LEGACY_ORPHAN_AGE_SECONDS
        if stale:
            shutil.rmtree(entry, ignore_errors=True)
            removed += 1
    return removed


def candidate_node_tmp_roots(
    environ: dict[str, str] | None = None,
    home_dir: Path | None = None,
) -> list[Path]:
    """Return writable node-local temp directories strictly outside $HOME."""
    env = os.environ if environ is None else environ
    home_resolved = str((home_dir or Path.home()).resolve())
    seen: set[str] = set()
    roots: list[Path] = []
    for raw in (env.get("SLURM_TMPDIR"), "/tmp", env.get("TMPDIR")):
        if not raw:
            continue
        p = Path(raw)
        try:
            if not p.is_dir() or not os.access(p, os.W_OK | os.X_OK):
                continue
            res = str(p.resolve())
        except OSError:
            continue
        if res == home_resolved or res.startswith(home_resolved + os.sep):
            continue
        if res not in seen:
            seen.add(res)
            roots.append(Path(res))
    return roots


def allocate_locus_workdir(
    safe_region: str,
    final_loci_root: Path,
    *,
    reserve_bytes: int = DEFAULT_LOCUS_RESERVE_BYTES,
    min_headroom_bytes: int = DEFAULT_MIN_HEADROOM_BYTES,
    tmp_candidates: list[Path] | None = None,
) -> WorkdirClaim:
    """Atomically reserve space on node-local /tmp or spill cleanly to `final_loci_root`."""
    uid = os.getuid()
    candidates = candidate_node_tmp_roots() if tmp_candidates is None else tmp_candidates

    for cand in candidates:
        lock_path = cand / f".gwasfm_space_{uid}.lock"
        state_path = cand / f".gwasfm_space_{uid}.json"
        try:
            with open(lock_path, "a+", encoding="utf-8") as lock_fh:
                fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
                _sweep_stale_tmp_dirs(cand, uid)

                raw_claims: dict[str, dict[str, Any]] = {}
                if state_path.is_file():
                    try:
                        loaded = json.loads(state_path.read_text(encoding="utf-8"))
                        if isinstance(loaded, dict) and isinstance(loaded.get("claims"), dict):
                            raw_claims = loaded["claims"]
                    except (OSError, ValueError, TypeError):
                        raw_claims = {}

                live_claims: dict[str, dict[str, Any]] = {}
                unwritten_reserved = 0
                for cid, info in raw_claims.items():
                    if not isinstance(info, dict):
                        continue
                    pid = int(info.get("pid", -1))
                    dpath = Path(str(info.get("dir", "")))
                    res_b = int(info.get("reserve", 0))
                    if _pid_is_alive(pid) and dpath.is_dir():
                        live_claims[cid] = {"pid": pid, "dir": str(dpath), "reserve": res_b}
                        written_b = _dir_size_bytes(dpath)
                        unwritten_reserved += max(0, res_b - written_b)

                free_now = shutil.disk_usage(cand).free
                effective_after = free_now - unwritten_reserved - reserve_bytes
                if effective_after >= min_headroom_bytes:
                    work_dir = Path(tempfile.mkdtemp(prefix=f"gwasfm_{safe_region}_", dir=str(cand)))
                    (work_dir / ".owner_pid").write_text(f"{os.getpid()}\n", encoding="utf-8")
                    claim_id = f"{os.getpid()}_{uuid.uuid4().hex[:8]}"
                    live_claims[claim_id] = {
                        "pid": os.getpid(),
                        "dir": str(work_dir),
                        "reserve": reserve_bytes,
                    }
                    tmp_state = state_path.with_suffix(".json.tmp")
                    tmp_state.write_text(json.dumps({"claims": live_claims}), encoding="utf-8")
                    os.replace(tmp_state, state_path)
                    return WorkdirClaim(
                        work_dir=work_dir,
                        using_node_tmp=True,
                        cleanup_on_exit=True,
                        lock_path=lock_path,
                        state_path=state_path,
                        claim_id=claim_id,
                    )
                else:
                    # Persist pruned live_claims even when spilling
                    tmp_state = state_path.with_suffix(".json.tmp")
                    tmp_state.write_text(json.dumps({"claims": live_claims}), encoding="utf-8")
                    os.replace(tmp_state, state_path)
        except OSError:
            continue

    final_loci_root.mkdir(parents=True, exist_ok=True)
    spill_dir = Path(tempfile.mkdtemp(prefix=f"gwasfm_{safe_region}_", dir=str(final_loci_root)))
    (spill_dir / ".owner_pid").write_text(f"{os.getpid()}\n", encoding="utf-8")
    return WorkdirClaim(
        work_dir=spill_dir,
        using_node_tmp=False,
        cleanup_on_exit=True,
    )


def allocate_explicit_temp_workdir(safe_region: str, explicit_temp_dir: Path) -> WorkdirClaim:
    """Stage locus work directly inside a caller-specified temp directory (e.g. /scratch)."""
    explicit_temp_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix=f"gwasfm_{safe_region}_", dir=str(explicit_temp_dir)))
    (work_dir / ".owner_pid").write_text(f"{os.getpid()}\n", encoding="utf-8")
    return WorkdirClaim(
        work_dir=work_dir,
        using_node_tmp=False,
        cleanup_on_exit=True,
    )


def release_locus_workdir(claim: WorkdirClaim) -> None:
    """Remove transient locus workdir and release any node-local /tmp space reservation."""
    if claim.cleanup_on_exit:
        shutil.rmtree(claim.work_dir, ignore_errors=True)
    if claim.using_node_tmp and claim.lock_path and claim.state_path and claim.claim_id:
        try:
            with open(claim.lock_path, "a+", encoding="utf-8") as lock_fh:
                fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
                if claim.state_path.is_file():
                    loaded = json.loads(claim.state_path.read_text(encoding="utf-8"))
                    claims = loaded.get("claims", {}) if isinstance(loaded, dict) else {}
                    if claim.claim_id in claims:
                        claims.pop(claim.claim_id, None)
                        tmp_state = claim.state_path.with_suffix(".json.tmp")
                        tmp_state.write_text(json.dumps({"claims": claims}), encoding="utf-8")
                        os.replace(tmp_state, claim.state_path)
        except (OSError, ValueError):
            pass


def is_enospc_failure(output_text: str) -> bool:
    low = (output_text or "").lower()
    return any(marker in low for marker in ENOSPC_MARKERS)


class StepExecutionError(RuntimeError):
    def __init__(self, returncode: int, cmd: list[str], combined_output: str) -> None:
        super().__init__(f"Command exited with status {returncode}: {' '.join(cmd)}")
        self.returncode = returncode
        self.cmd = cmd
        self.combined_output = combined_output


def _run_and_log(
    cmd: list[str],
    step_stem: str,
    gwas_log_dir: Path,
    sub_env: dict[str, str],
) -> None:
    res = subprocess.run(cmd, text=True, capture_output=True, env=sub_env)
    out_txt = res.stdout or ""
    err_txt = res.stderr or ""
    combined = out_txt + (("\n" + err_txt) if err_txt and res.returncode != 0 else "")
    log_path = gwas_log_dir / f"{step_stem}.log"
    try:
        log_path.write_text(combined if combined.strip() else f"[completed] {step_stem}\n", encoding="utf-8")
        if err_txt.strip():
            (gwas_log_dir / f"{step_stem}.err").write_text(err_txt, encoding="utf-8")
    except OSError as exc:
        sys.stderr.write(f"[WARN] Could not write step log {log_path}: {exc}\n")
    if res.returncode != 0:
        raise StepExecutionError(res.returncode, cmd, combined)


def _build_prep_cmd(
    proj: Path,
    row: dict[str, Any],
    screen_json_path: str,
    loc_idx: str,
    work_loc_dir: Path,
    keep_inter: bool,
) -> list[str]:
    prep_cache = str(row["cache_dir"]) if keep_inter else str(work_loc_dir / "cache")
    cmd = [
        sys.executable,
        str(proj / "scripts" / "prep_locus.py"),
        "--screen-json",
        screen_json_path,
        "--locus-index",
        str(loc_idx),
        "--outdir",
        str(work_loc_dir),
        "--chain",
        row["chain"],
        "--grch37-fasta",
        row["grch37_fasta"],
        "--cache-dir",
        prep_cache,
    ]
    if row.get("panukb_dir"):
        cmd.extend(["--panukb-dir", row["panukb_dir"]])
    if row.get("local_ld_panel"):
        cmd.extend(["--local-ld-panel", row["local_ld_panel"]])
    if row.get("bm_path"):
        cmd.extend(["--bm-path", row["bm_path"]])
    if keep_inter:
        cmd.append("--keep-intermediates")
    return cmd


def _make_sub_env(work_loc_dir: Path) -> dict[str, str]:
    work_tmp_dir = work_loc_dir / "tmp"
    work_tmp_dir.mkdir(parents=True, exist_ok=True)
    sub_env = dict(os.environ)
    sub_env["TMPDIR"] = str(work_tmp_dir)
    return sub_env


def run_locus_worker(args: argparse.Namespace) -> int:
    proj = Path(args.project_dir).resolve()
    row = json.loads(Path(args.row_json).read_text(encoding="utf-8"))
    screen = json.loads(Path(args.screen_json).read_text(encoding="utf-8"))
    loc_idx = str(args.locus_index)
    reg_id = args.region_id
    reg_type = args.region_type
    scratch_base = Path(args.scratch_base)

    keep_inter = bool(row.get("keep_intermediates"))
    gwas_log_dir = Path(row["outdir"]) / row["study_id"] / row["gwas_id"] / "logs"
    gwas_log_dir.mkdir(parents=True, exist_ok=True)
    safe_reg = re.sub(r"[^A-Za-z0-9._-]+", "_", reg_id)

    if keep_inter:
        final_loci_root = Path(row["outdir"]) / row["study_id"] / row["gwas_id"] / "loci"
        final_loc_dir = final_loci_root / reg_id
        final_loc_dir.mkdir(parents=True, exist_ok=True)
        claim = WorkdirClaim(
            work_dir=final_loc_dir,
            using_node_tmp=False,
            cleanup_on_exit=False,
        )
    else:
        final_loci_root = scratch_base
        final_loc_dir = final_loci_root / reg_id
        final_loc_dir.mkdir(parents=True, exist_ok=True)
        explicit_temp = (args.temp_dir or "").strip()
        if explicit_temp:
            claim = allocate_explicit_temp_workdir(safe_reg, Path(explicit_temp).expanduser())
        else:
            reserve_b = MHC_LOCUS_RESERVE_BYTES if reg_type == "mhc" else DEFAULT_LOCUS_RESERVE_BYTES
            claim = allocate_locus_workdir(safe_reg, final_loci_root, reserve_bytes=reserve_b)

    try:
        sub_env = _make_sub_env(claim.work_dir)
        prep_cmd = _build_prep_cmd(proj, row, args.screen_json, loc_idx, claim.work_dir, keep_inter)
        try:
            _run_and_log(prep_cmd, f"02_prep_{safe_reg}", gwas_log_dir, sub_env)
        except StepExecutionError as exc:
            if claim.using_node_tmp and is_enospc_failure(exc.combined_output):
                sys.stderr.write(
                    f"[WARN] Node-local temp filled during {reg_id}; spilling locus prep to {final_loci_root}\n"
                )
                release_locus_workdir(claim)
                claim = allocate_locus_workdir(
                    safe_reg,
                    final_loci_root,
                    tmp_candidates=[],
                )
                sub_env = _make_sub_env(claim.work_dir)
                prep_cmd = _build_prep_cmd(proj, row, args.screen_json, loc_idx, claim.work_dir, keep_inter)
                _run_and_log(prep_cmd, f"02_prep_{safe_reg}", gwas_log_dir, sub_env)
            else:
                sys.stderr.write(exc.combined_output)
                if exc.returncode < 0:
                    sig = abs(exc.returncode)
                    sys.stderr.write(f"[ERROR] Locus prep killed by signal {sig} (likely OOM): {' '.join(exc.cmd)}\n")
                    return 128 + sig
                raise subprocess.CalledProcessError(exc.returncode, exc.cmd) from exc

        prep_qc_src = claim.work_dir / "locus_prep_qc.json"
        if prep_qc_src.is_file() and claim.work_dir != final_loc_dir:
            shutil.copy2(prep_qc_src, final_loc_dir / "locus_prep_qc.json")

        prep_qc = json.loads(prep_qc_src.read_text(encoding="utf-8")) if prep_qc_src.is_file() else {}
        if prep_qc.get("status") == "prepared":
            variants_tsv = claim.work_dir / "matched_variants.tsv"
            matrix_npy = claim.work_dir / "matrix.npy"
            for method_id in row["methods"]:
                script_path = proj / "scripts" / "methods" / f"run_{method_id}.R"
                if not script_path.is_file():
                    raise FileNotFoundError(f"Standalone method script not found: {script_path}")
                out_m_dir = final_loc_dir / "methods" / method_id
                out_m_dir.mkdir(parents=True, exist_ok=True)
                m_cmd = [
                    "Rscript",
                    str(script_path),
                    "--variants",
                    str(variants_tsv),
                    "--matrix",
                    str(matrix_npy),
                    "--region-id",
                    reg_id,
                    "--region-type",
                    reg_type,
                    "--out-rds",
                    str(out_m_dir / "region_result.rds"),
                    "--qc-json",
                    str(out_m_dir / "method_qc.json"),
                    "--allow-skipped-locus",
                    "true",
                    "--compress",
                    "false" if keep_inter else "gzip",
                ]
                if method_id != "p_ld":
                    m_cmd.extend([
                        "--max-causal",
                        str(row.get("max_causal", 10)),
                        "--max-iter",
                        str(row.get("max_iter", 50)),
                        "--reference-sample-size",
                        str(screen["reference_sample_size"]),
                    ])
                    if row.get("effective_n") is not None:
                        m_cmd.extend(["--effective-n", str(row["effective_n"])])
                    if method_id == "polyfun" and row.get("polyfun_dir"):
                        m_cmd.extend(["--polyfun-priors", str(row["polyfun_dir"])])

                try:
                    _run_and_log(m_cmd, f"03_method_{safe_reg}_{method_id}", gwas_log_dir, sub_env)
                except StepExecutionError as exc:
                    sys.stderr.write(exc.combined_output)
                    if exc.returncode < 0:
                        sig = abs(exc.returncode)
                        sys.stderr.write(f"[ERROR] Method {method_id} killed by signal {sig} (likely OOM)\n")
                        return 128 + sig
                    raise subprocess.CalledProcessError(exc.returncode, exc.cmd) from exc
    finally:
        release_locus_workdir(claim)

    Path(args.out_root_file).write_text(str(final_loci_root), encoding="utf-8")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--row-json", required=True)
    p.add_argument("--screen-json", required=True)
    p.add_argument("--locus-index", required=True)
    p.add_argument("--region-id", required=True)
    p.add_argument("--region-type", required=True)
    p.add_argument("--scratch-base", required=True)
    p.add_argument("--temp-dir", default="")
    p.add_argument("--out-root-file", default=".scratch_loci_root")
    return p.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(run_locus_worker(parse_args()))
