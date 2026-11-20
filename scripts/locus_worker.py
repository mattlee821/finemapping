#!/usr/bin/env python3
"""Per-locus worker for Nextflow FINEMAP_LOCUS with per-locus checkpointing and zero-NFS staging.

1. Keys completed per-locus `.rds` / `method_qc.json` outputs deterministically under
   `<scratch_root>/by_gwas/<study_id>__<gwas_id>__<settings_hash>/loci/<region_id>/`
   so interrupted runs immediately reuse already-completed loci on restart.
2. Coordinates transient multi-GB LD matrix (`matrix.npy`) staging across concurrent
   workers on the same compute node via an `flock`-guarded reservation table in
   `$SLURM_TMPDIR` -> `/dev/shm` (RAM `tmpfs`) -> `/usr/local/scratch` -> `/tmp`
   using dynamic $O(M^2)$ byte reservations, keeping multi-GB `matrix.npy` writes off
   shared `/scratch` NFS (`balance_dirty_pages`).
3. Sweeps orphaned `gwasfm_*` temp directories from cancelled jobs on the local host,
   runs `prep_locus.py` followed by `scripts/methods/run_all_methods.R` (which unlinks
   `matrix.npy` immediately after reading), and releases the workdir in `finally:`.
"""

from __future__ import annotations

import argparse
import errno
import fcntl
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_LOCUS_RESERVE_BYTES = int(5.5 * (1024 ** 3))   # Fallback for standard <=9.5 Mb locus
MHC_LOCUS_RESERVE_BYTES = int(55.0 * (1024 ** 3))      # Fallback for extended 8.5 Mb MHC locus (~80k variants)
MIN_DYNAMIC_RESERVE_BYTES = int(1.5 * (1024 ** 3))     # Floor for sparse loci
DEFAULT_MIN_HEADROOM_BYTES = int(4.0 * (1024 ** 3))    # Leave >=4 GiB untouched on node temp/tmpfs
LEGACY_ORPHAN_AGE_SECONDS = 120.0
SKIPPED_PREP_STATUSES = ("skipped_no_lifted_variants", "skipped_no_panel_matched_variants")

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
    """Delete orphaned gwasfm_* directories on THIS host owned by `uid` whose creator PID exited."""
    now_ts = time.time() if now is None else now
    my_host = socket.gethostname()
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
                raw_owner = pid_file.read_text(encoding="utf-8").strip()
                if ":" in raw_owner:
                    owner_host, pid_str = raw_owner.rsplit(":", 1)
                    if owner_host != my_host:
                        # Created by a worker on a different compute node — never delete!
                        continue
                    stale = not _pid_is_alive(int(pid_str))
                else:
                    owner_pid = int(raw_owner)
                    stale = not _pid_is_alive(owner_pid)
            except (OSError, ValueError):
                stale = False
        else:
            stale = (now_ts - st.st_mtime) >= 3600.0
        if stale:
            shutil.rmtree(entry, ignore_errors=True)
            removed += 1
    return removed


def candidate_node_tmp_roots(
    environ: dict[str, str] | None = None,
    home_dir: Path | None = None,
    exclude_roots: list[Path] | None = None,
) -> list[Path]:
    """Return writable node-local temp directories strictly outside $HOME and shared /scratch."""
    env = os.environ if environ is None else environ
    home_resolved = str((home_dir or Path.home()).resolve())
    excluded_prefixes = [home_resolved, "/scratch"]
    for ex in (exclude_roots or []):
        try:
            excluded_prefixes.append(str(ex.resolve()))
        except OSError:
            pass

    seen: set[str] = set()
    roots: list[Path] = []
    uid = os.getuid()
    user_name = env.get("USER") or str(uid)
    # Order: dedicated SLURM job SSD -> RAM tmpfs (/dev/shm) -> node-local disk (/usr/local/scratch) -> /tmp.
    # Never include $TMPDIR here because finemap.sh exports TMPDIR to shared /scratch/.../tmp/<RUN_ID>.
    raw_candidates = (
        env.get("SLURM_TMPDIR"),
        "/dev/shm",
        f"/usr/local/scratch/{user_name}",
        "/usr/local/scratch",
        "/tmp",
    )
    for raw in raw_candidates:
        if not raw:
            continue
        p = Path(raw)
        try:
            if not p.is_dir() or not os.access(p, os.W_OK | os.X_OK):
                continue
            res = str(p.resolve())
        except OSError:
            continue
        if any(res == pref or res.startswith(pref.rstrip(os.sep) + os.sep) for pref in excluded_prefixes if pref):
            continue
        if res not in seen:
            seen.add(res)
            roots.append(Path(res))
    return roots


def estimate_locus_reserve_bytes(screen: dict[str, Any], loc_idx: int, reg_type: str) -> int:
    """Estimate transient workdir bytes needed for `matrix.npy` + lifted VCF/index from variant count M."""
    regions = screen.get("selected_regions") or []
    v_cnt: int | None = None
    if 0 <= loc_idx < len(regions) and isinstance(regions[loc_idx], dict):
        raw_v = regions[loc_idx].get("vcf_variants_in_region")
        if raw_v is not None:
            try:
                v_cnt = max(1, int(raw_v))
            except (TypeError, ValueError):
                v_cnt = None
    if v_cnt is None:
        return MHC_LOCUS_RESERVE_BYTES if reg_type == "mhc" else DEFAULT_LOCUS_RESERVE_BYTES
    matrix_bytes = v_cnt * v_cnt * 8
    return max(MIN_DYNAMIC_RESERVE_BYTES, int(matrix_bytes * 1.15) + 512 * (1024 ** 2))


def gwas_scratch_loci_dir(scratch_root: Path, row: dict[str, Any]) -> Path:
    """Deterministic per-(study_id, gwas_id, settings_hash) scratch root for completed loci."""
    study_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(row.get("study_id") or "study"))[:80]
    gwas_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(row.get("gwas_id") or "gwas"))[:80]
    shash = re.sub(r"[^A-Za-z0-9._-]+", "_", str(row.get("settings_hash") or "default"))[:32]
    return scratch_root / "by_gwas" / f"{study_safe}__{gwas_safe}__{shash}" / "loci"


def _normalize_method_id(name: str) -> str:
    m = name.strip().lower().replace("-", "_")
    return "polyfun" if m == "polyfun_susie" else m


def locus_outputs_are_complete(loc_dir: Path, methods: list[str]) -> bool:
    """Return True iff `loc_dir` already holds valid prep QC and every requested method output."""
    prep_qc_path = loc_dir / "locus_prep_qc.json"
    try:
        if not prep_qc_path.is_file() or prep_qc_path.stat().st_size == 0:
            return False
        prep_qc = json.loads(prep_qc_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return False
    if not isinstance(prep_qc, dict):
        return False

    status = prep_qc.get("status")
    if status in SKIPPED_PREP_STATUSES:
        return True
    if status != "prepared":
        return False

    for m_raw in methods:
        m_id = _normalize_method_id(m_raw)
        rds_path = loc_dir / "methods" / m_id / "region_result.rds"
        mqc_path = loc_dir / "methods" / m_id / "method_qc.json"
        try:
            if not rds_path.is_file() or rds_path.stat().st_size == 0:
                return False
            if not mqc_path.is_file() or mqc_path.stat().st_size == 0:
                return False
            mqc = json.loads(mqc_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            return False
        if not isinstance(mqc, dict) or mqc.get("method_id") != m_id:
            return False
    return True


def migrate_run_scratch(run_scratch_dir: Path, rows_dir: Path, scratch_root: Path) -> dict[str, int]:
    """Migrate completed locus outputs from a legacy `scratch/<RUN_ID>/row_*/loci/*` tree into `by_gwas/`."""
    migrated_loci = 0
    touched_gwas = 0
    if not run_scratch_dir.is_dir() or not rows_dir.is_dir():
        return {"migrated_loci": 0, "gwas_with_migrated_loci": 0}

    for row_json in sorted(rows_dir.glob("row_*.json")):
        try:
            row = json.loads(row_json.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        old_loci = run_scratch_dir / row_json.stem / "loci"
        if not old_loci.is_dir():
            continue
        methods = row.get("methods") or ["p_ld", "susie_rss", "susie_rrss", "carma", "polyfun"]
        dst_loci = gwas_scratch_loci_dir(scratch_root, row)
        gwas_moved = 0
        for loc_dir in sorted(old_loci.iterdir()):
            if not loc_dir.is_dir() or loc_dir.name.startswith("gwasfm_"):
                continue
            if not locus_outputs_are_complete(loc_dir, methods):
                continue
            dst_loc = dst_loci / loc_dir.name
            if locus_outputs_are_complete(dst_loc, methods):
                continue
            dst_loci.mkdir(parents=True, exist_ok=True)
            if dst_loc.exists():
                shutil.rmtree(dst_loc, ignore_errors=True)
            shutil.move(str(loc_dir), str(dst_loc))
            migrated_loci += 1
            gwas_moved += 1
        if gwas_moved > 0:
            touched_gwas += 1
    return {"migrated_loci": migrated_loci, "gwas_with_migrated_loci": touched_gwas}


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
    my_host = socket.gethostname()
    candidates = (
        candidate_node_tmp_roots(exclude_roots=[final_loci_root])
        if tmp_candidates is None
        else tmp_candidates
    )

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
                    host = str(info.get("host", my_host))
                    dpath = Path(str(info.get("dir", "")))
                    res_b = int(info.get("reserve", 0))
                    alive = (host != my_host) or _pid_is_alive(pid)
                    if alive and dpath.is_dir():
                        live_claims[cid] = {"host": host, "pid": pid, "dir": str(dpath), "reserve": res_b}
                        written_b = _dir_size_bytes(dpath)
                        unwritten_reserved += max(0, res_b - written_b)

                free_now = shutil.disk_usage(cand).free
                effective_after = free_now - unwritten_reserved - reserve_bytes
                if effective_after >= min_headroom_bytes:
                    work_dir = Path(tempfile.mkdtemp(prefix=f"gwasfm_{safe_region}_", dir=str(cand)))
                    (work_dir / ".owner_pid").write_text(f"{my_host}:{os.getpid()}\n", encoding="utf-8")
                    claim_id = f"{my_host}_{os.getpid()}_{uuid.uuid4().hex[:8]}"
                    live_claims[claim_id] = {
                        "host": my_host,
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
    (spill_dir / ".owner_pid").write_text(f"{my_host}:{os.getpid()}\n", encoding="utf-8")
    return WorkdirClaim(
        work_dir=spill_dir,
        using_node_tmp=False,
        cleanup_on_exit=True,
    )


def allocate_explicit_temp_workdir(safe_region: str, explicit_temp_dir: Path) -> WorkdirClaim:
    """Stage locus work directly inside a caller-specified temp directory (e.g. /scratch)."""
    explicit_temp_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(tempfile.mkdtemp(prefix=f"gwasfm_{safe_region}_", dir=str(explicit_temp_dir)))
    (work_dir / ".owner_pid").write_text(f"{socket.gethostname()}:{os.getpid()}\n", encoding="utf-8")
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
    bm_disk_cache: str = "",
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
    if bm_disk_cache:
        cmd.extend(["--bm-disk-cache", bm_disk_cache])
    if keep_inter:
        cmd.append("--keep-intermediates")
    return cmd


def _make_sub_env(work_loc_dir: Path) -> dict[str, str]:
    work_tmp_dir = work_loc_dir / "tmp"
    work_tmp_dir.mkdir(parents=True, exist_ok=True)
    sub_env = dict(os.environ)
    sub_env["TMPDIR"] = str(work_tmp_dir)
    return sub_env


def _record_final_locus_failure(
    final_loc_dir: Path,
    reg_id: str,
    reg_type: str,
    error_msg: str,
    out_root_file: str,
    final_loci_root: Path,
) -> int:
    final_loc_dir.mkdir(parents=True, exist_ok=True)
    fail_qc = {
        "region_id": reg_id,
        "region_type": reg_type,
        "status": "failed_locus",
        "error": error_msg[:1000],
    }
    (final_loc_dir / "locus_prep_qc.json").write_text(
        json.dumps(fail_qc, indent=2) + "\n",
        encoding="utf-8",
    )
    Path(out_root_file).write_text(str(final_loci_root), encoding="utf-8")
    return 0


def run_locus_worker(args: argparse.Namespace) -> int:
    if args.migrate_run_scratch:
        if not args.rows_dir or not args.scratch_root:
            raise ValueError("--migrate-run-scratch requires both --rows-dir and --scratch-root")
        stats = migrate_run_scratch(
            Path(args.migrate_run_scratch).expanduser().resolve(),
            Path(args.rows_dir).expanduser().resolve(),
            Path(args.scratch_root).expanduser().resolve(),
        )
        print(
            f"[migrate_run_scratch] Migrated {stats['migrated_loci']} completed loci "
            f"across {stats['gwas_with_migrated_loci']} GWAS into {Path(args.scratch_root) / 'by_gwas'}"
        )
        return 0

    proj = Path(args.project_dir).resolve()
    row = json.loads(Path(args.row_json).read_text(encoding="utf-8"))
    screen = json.loads(Path(args.screen_json).read_text(encoding="utf-8"))
    loc_idx = str(args.locus_index)
    loc_idx_int = int(args.locus_index)
    reg_id = args.region_id
    reg_type = args.region_type
    bm_disk_cache = (args.bm_disk_cache or os.environ.get("GWASFM_BM_CACHE_DIR", "")).strip()
    is_final_attempt = int(args.attempt) >= int(args.max_attempts)

    keep_inter = bool(row.get("keep_intermediates"))
    force_run = bool(row.get("force"))
    gwas_log_dir = Path(row["outdir"]) / row["study_id"] / row["gwas_id"] / "logs"
    gwas_log_dir.mkdir(parents=True, exist_ok=True)
    safe_reg = re.sub(r"[^A-Za-z0-9._-]+", "_", reg_id)

    if keep_inter:
        final_loci_root = Path(row["outdir"]) / row["study_id"] / row["gwas_id"] / "loci"
    elif args.scratch_root:
        final_loci_root = gwas_scratch_loci_dir(Path(args.scratch_root), row)
    elif args.scratch_base:
        final_loci_root = Path(args.scratch_base)
    else:
        raise ValueError("Supply either --scratch-root or --scratch-base")

    final_loc_dir = final_loci_root / reg_id

    # Fast per-locus checkpoint skip when complete outputs already exist from a prior run
    if (not keep_inter) and (not force_run) and locus_outputs_are_complete(final_loc_dir, row["methods"]):
        try:
            (gwas_log_dir / f"02_prep_{safe_reg}.reuse.log").write_text(
                f"[locus_worker] Reusing completed locus {reg_id} from {final_loc_dir}\n",
                encoding="utf-8",
            )
        except OSError:
            pass
        Path(args.out_root_file).write_text(str(final_loci_root), encoding="utf-8")
        return 0

    final_loc_dir.mkdir(parents=True, exist_ok=True)
    if not keep_inter:
        # Clear any partial/failed state for THIS locus before re-running
        (final_loc_dir / "locus_prep_qc.json").unlink(missing_ok=True)
        shutil.rmtree(final_loc_dir / "methods", ignore_errors=True)
    for stale_log in gwas_log_dir.glob(f"0[23]_*_{safe_reg}*"):
        stale_log.unlink(missing_ok=True)

    if keep_inter:
        spill_target = final_loci_root
        claim = WorkdirClaim(
            work_dir=final_loc_dir,
            using_node_tmp=False,
            cleanup_on_exit=False,
        )
    else:
        explicit_temp = (args.temp_dir or "").strip()
        spill_target = Path(explicit_temp).expanduser() if explicit_temp else final_loci_root
        reserve_b = estimate_locus_reserve_bytes(screen, loc_idx_int, reg_type)
        # Stage in node-local $SLURM_TMPDIR / /dev/shm / /usr/local/scratch / /tmp first (zero NFS traffic)
        claim = allocate_locus_workdir(safe_reg, spill_target, reserve_bytes=reserve_b)

    try:
        sub_env = _make_sub_env(claim.work_dir)
        prep_cmd = _build_prep_cmd(proj, row, args.screen_json, loc_idx, claim.work_dir, keep_inter, bm_disk_cache)
        try:
            _run_and_log(prep_cmd, f"02_prep_{safe_reg}", gwas_log_dir, sub_env)
        except StepExecutionError as exc:
            if claim.using_node_tmp and is_enospc_failure(exc.combined_output):
                sys.stderr.write(
                    f"[WARN] Node-local temp filled during {reg_id}; spilling locus prep to {spill_target}\n"
                )
                release_locus_workdir(claim)
                claim = allocate_locus_workdir(
                    safe_reg,
                    spill_target,
                    tmp_candidates=[],
                )
                sub_env = _make_sub_env(claim.work_dir)
                prep_cmd = _build_prep_cmd(
                    proj, row, args.screen_json, loc_idx, claim.work_dir, keep_inter, bm_disk_cache
                )
                try:
                    _run_and_log(prep_cmd, f"02_prep_{safe_reg}", gwas_log_dir, sub_env)
                except StepExecutionError as exc2:
                    sys.stderr.write(exc2.combined_output)
                    if is_final_attempt:
                        return _record_final_locus_failure(
                            final_loc_dir, reg_id, reg_type, f"prep_locus failed: {exc2}", args.out_root_file, final_loci_root
                        )
                    if exc2.returncode < 0:
                        return 128 + abs(exc2.returncode)
                    raise subprocess.CalledProcessError(exc2.returncode, exc2.cmd) from exc2
            else:
                sys.stderr.write(exc.combined_output)
                if is_final_attempt:
                    return _record_final_locus_failure(
                        final_loc_dir, reg_id, reg_type, f"prep_locus failed: {exc}", args.out_root_file, final_loci_root
                    )
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
            batch_runner = proj / "scripts" / "methods" / "run_all_methods.R"
            if not batch_runner.is_file():
                raise FileNotFoundError(f"Batch method runner not found: {batch_runner}")
            m_cmd = [
                "Rscript",
                str(batch_runner),
                "--variants",
                str(variants_tsv),
                "--matrix",
                str(matrix_npy),
                "--region-id",
                reg_id,
                "--region-type",
                reg_type,
                "--methods",
                ",".join(row["methods"]),
                "--out-methods-dir",
                str(final_loc_dir / "methods"),
                "--log-dir",
                str(gwas_log_dir),
                "--safe-region",
                safe_reg,
                "--max-causal",
                str(row.get("max_causal", 10)),
                "--max-iter",
                str(row.get("max_iter", 50)),
                "--reference-sample-size",
                str(screen["reference_sample_size"]),
                "--allow-skipped-locus",
                "true",
                "--compress",
                "false" if keep_inter else "gzip",
                "--unlink-matrix",
                "false" if keep_inter else "true",
            ]
            if row.get("effective_n") is not None:
                m_cmd.extend(["--effective-n", str(row["effective_n"])])
            if row.get("polyfun_dir"):
                m_cmd.extend(["--polyfun-priors", str(row["polyfun_dir"])])

            try:
                _run_and_log(m_cmd, f"03_methods_batch_{safe_reg}", gwas_log_dir, sub_env)
            except StepExecutionError as exc:
                sys.stderr.write(exc.combined_output)
                if is_final_attempt:
                    return _record_final_locus_failure(
                        final_loc_dir, reg_id, reg_type, f"methods failed: {exc}", args.out_root_file, final_loci_root
                    )
                if exc.returncode < 0:
                    sig = abs(exc.returncode)
                    sys.stderr.write(f"[ERROR] Methods runner killed by signal {sig} (likely OOM)\n")
                    return 128 + sig
                raise subprocess.CalledProcessError(exc.returncode, exc.cmd) from exc
    finally:
        release_locus_workdir(claim)

    Path(args.out_root_file).write_text(str(final_loci_root), encoding="utf-8")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-dir", default="")
    p.add_argument("--row-json", default="")
    p.add_argument("--screen-json", default="")
    p.add_argument("--locus-index", default="0")
    p.add_argument("--region-id", default="")
    p.add_argument("--region-type", default="lead_locus")
    p.add_argument("--scratch-root", default="")
    p.add_argument("--scratch-base", default="")
    p.add_argument("--temp-dir", default="")
    p.add_argument("--bm-disk-cache", default="")
    p.add_argument("--attempt", type=int, default=1)
    p.add_argument("--max-attempts", type=int, default=1)
    p.add_argument("--out-root-file", default=".scratch_loci_root")
    p.add_argument("--migrate-run-scratch", default="", help="Migrate completed loci from scratch/<RUN_ID> into by_gwas/")
    p.add_argument("--rows-dir", default="", help="Rows JSON directory matching --migrate-run-scratch")
    return p.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(run_locus_worker(parse_args()))
