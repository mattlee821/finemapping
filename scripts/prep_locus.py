#!/usr/bin/env python3
"""
Standalone Stage 2: Shared per-locus data preparation (liftover, panel match, LD extraction).

Produces shared inputs consumed by all standalone fine-mapping method scripts:
  - <locus_dir>/matched_variants.tsv   (aligned GRCh38 + GRCh37 table with EZ and LD_R2_TO_LEAD)
  - <locus_dir>/matrix.npy             (Fortran-order float64 normalized symmetric LD matrix)
  - <locus_dir>/matrix.ids.txt         (1-line-per-variant companion IDs)
  - <locus_dir>/locus_prep_qc.json     (variant attrition, LD validation, SS stats, sub-step timings)

Can be run:
  1. Interactively / with no arguments (uses the default local NPPB cis_NPPB locus & local test BM if present), OR
  2. From CLI / Nextflow:
     python3 scripts/prep_locus.py \
       --screen-json /path/to/screen_summary.json \
       --locus-index 0 \
       --outdir /path/to/loci/cis_NPPB_chr1_10857464_12858945
"""

from __future__ import annotations

import argparse
from collections import Counter, OrderedDict
import gc
import json
import os
from pathlib import Path
import random
import resource
import sys
import tempfile
import threading
import time
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_finemapping import (  # noqa: E402
    MHC_CHROM,
    MHC_CORE_GRC_GRCH37,
    _normalize_ld_matrix_inplace,
    _write_symmetric_fortran_npy,
    chromosome_sort_key,
    default_panukb_dir,
    ensure_regional_index,
    extract_ld_from_local_panel,
    format_vcf_region,
    inspect_vcf_tags,
    is_mhc_core_grc_grch38,
    is_mhc_grch38,
    match_panel_variants,
    matrix_to_npy,
    normalize_chromosome,
    normalize_lifted_vcf,
    query_vcf,
    resolve_panel_matches,
    write_json,
    write_lifted_vcf,
    write_tsv,
)

# ==============================================================================
# 1. Interactive defaults (used when run without CLI arguments)
# ==============================================================================
HOME = Path.home()
DEFAULT_SCREEN_JSON = HOME / "research-local/GWASfinemapping/test_nppb/screen_summary.json"
DEFAULT_REF_ROOT = Path(os.environ.get("GWASFM_REFERENCE_ROOT", str(HOME / ".cache/gwasfinemapping/reference")))
DEFAULT_CACHE_DIR = Path(os.environ.get("GWASFM_CACHE_DIR", str(HOME / ".cache/gwasfinemapping")))
DEFAULT_LOCAL_TEST_BM = HOME / "research-local/GWASfinemapping/nppb_test_bm/UKBB.EUR.ldadj.bm"
DEFAULT_LOCUS_ROOT = HOME / "research-local/GWASfinemapping/test_nppb/loci"


def current_peak_rss_mb() -> float:
    ru_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ru_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak = max(ru_self, ru_child)
    if sys.platform == "darwin":
        return round(peak / (1024.0 * 1024.0), 1)
    return round(peak / 1024.0, 1)


def current_cpu_times() -> tuple[float, float]:
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return round(s.ru_utime + c.ru_utime, 3), round(s.ru_stime + c.ru_stime, 3)


def open_blockmatrix_reader(
    bm_uri: str,
    block_cache: int = 4,
    disk_cache_root: str | None = None,
) -> Any:
    """Open a Hail BlockMatrix reader with optional atomic local-disk part-* caching."""
    from ldcov.io.blockmatrix.block_codec import decode_block
    from ldcov.io.blockmatrix.partitioner import GridPartitioner
    from ldcov.io.blockmatrix.reader import HailBlockMatrixReader
    from ldcov.io.fs_utils import resolve_filesystem

    cache_str = (disk_cache_root or "").strip()
    is_remote = bm_uri.startswith(("s3://", "gs://", "http://", "https://"))
    if not cache_str or not is_remote:
        reader = HailBlockMatrixReader(bm_uri, block_cache=block_cache)
        reader.disk_hits = 0
        reader.disk_misses = 0
        return reader

    class DiskCachedHailBlockMatrixReader:
        def __init__(self, path: str, cache_root: Path, mem_cache_size: int = 4) -> None:
            self.path = path.rstrip("/")
            self._fs = None
            self._fs_lock = threading.Lock()
            self._cache_size = max(1, mem_cache_size)
            self._cache: OrderedDict[tuple[int, int], Any] = OrderedDict()
            self._lock = threading.Lock()
            self.disk_hits = 0
            self.disk_misses = 0

            bm_name = self.path.rsplit("/", 1)[-1]
            self.disk_dir = cache_root / bm_name
            meta_local = self.disk_dir / "metadata.json"
            meta = None
            if meta_local.is_file() and meta_local.stat().st_size > 10:
                try:
                    meta = json.loads(meta_local.read_text(encoding="utf-8"))
                except Exception:
                    meta = None
            if meta is None:
                raw_meta = self._fetch_remote_bytes(self.path + "/metadata.json")
                meta = json.loads(raw_meta.decode("utf-8"))
                self._atomic_cache_write(meta_local, raw_meta)

            self.block_size = meta["blockSize"]
            self.n_rows = meta["nRows"]
            self.n_cols = meta["nCols"]
            self._part_files = meta["partFiles"]
            self.partitioner = GridPartitioner(
                self.block_size, self.n_rows, self.n_cols, meta.get("maybeFiltered")
            )

        def _get_fs(self) -> Any:
            if self._fs is None:
                with self._fs_lock:
                    if self._fs is None:
                        self._fs, _ = resolve_filesystem(self.path, None)
            return self._fs

        def _fetch_remote_bytes(self, remote_path: str, max_attempts: int = 4) -> bytes:
            last_exc: Exception | None = None
            for attempt in range(max_attempts):
                try:
                    fs = self._get_fs()
                    with fs.open(remote_path, "rb") as fh:
                        return fh.read()
                except Exception as exc:  # noqa: BLE001
                    last_exc = exc
                    if attempt + 1 < max_attempts:
                        time.sleep(0.5 * (2 ** attempt) + random.uniform(0.0, 0.35))
                        with self._fs_lock:
                            self._fs = None
            assert last_exc is not None
            raise last_exc

        @staticmethod
        def _atomic_cache_write(dst: Path, payload: bytes) -> None:
            try:
                dst.parent.mkdir(parents=True, exist_ok=True)
                tmp = dst.parent / f".{dst.name}.tmp.{os.getpid()}.{threading.get_ident()}"
                with open(tmp, "wb") as fh:
                    fh.write(payload)
                    fh.flush()
                    try:
                        os.fsync(fh.fileno())
                    except OSError:
                        pass
                os.replace(tmp, dst)
            except OSError:
                pass

        def read_block(self, i: int, j: int) -> Any:
            key = (i, j)
            with self._lock:
                if key in self._cache:
                    self._cache.move_to_end(key)
                    return self._cache[key]

            slot = self.partitioner.part_slot(i, j)
            if slot is None:
                return None

            part_name = self._part_files[slot]
            local_part = self.disk_dir / "parts" / part_name
            block = None
            if local_part.is_file() and local_part.stat().st_size > 0:
                try:
                    block = decode_block(local_part.read_bytes())
                    with self._lock:
                        self.disk_hits += 1
                except Exception:  # noqa: BLE001
                    local_part.unlink(missing_ok=True)
                    block = None

            if block is None:
                raw = self._fetch_remote_bytes(self.path + "/parts/" + part_name)
                block = decode_block(raw)
                self._atomic_cache_write(local_part, raw)
                with self._lock:
                    self.disk_misses += 1

            with self._lock:
                self._cache[key] = block
                self._cache.move_to_end(key)
                if len(self._cache) > self._cache_size:
                    self._cache.popitem(last=False)
            return block

    return DiskCachedHailBlockMatrixReader(
        bm_uri,
        Path(cache_str).expanduser().resolve(),
        mem_cache_size=block_cache,
    )


def _assemble_bounded(reader: Any, global_idxs: list[int], max_workers: int = 4) -> Any:
    """Assemble symmetric LD submatrix with bounded in-flight blocks (prevents O(blocks) RAM spike on MHC)."""
    from collections import defaultdict
    from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
    import numpy as np
    from ldcov.ld_bm.extract import _as_slice

    bs = reader.block_size
    n = len(global_idxs)
    out = np.full((n, n), np.nan, dtype=np.float64)

    groups: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for pos, g in enumerate(global_idxs):
        groups[g // bs].append((pos, g % bs))
    block_ids = sorted(groups)

    pairs = [
        (p, q)
        for a_i, p in enumerate(block_ids)
        for q in block_ids[a_i:]
        if reader.partitioner.part_slot(p, q) is not None
    ]

    def _scatter(p: int, q: int, block: Any) -> None:
        a_positions = [pos for pos, _ in groups[p]]
        a_locals = [loc for _, loc in groups[p]]
        b_positions = [pos for pos, _ in groups[q]]
        b_locals = [loc for _, loc in groups[q]]
        a_ps, a_ls = _as_slice(a_positions), _as_slice(a_locals)
        b_ps, b_ls = _as_slice(b_positions), _as_slice(b_locals)
        if a_ps and a_ls and b_ps and b_ls:
            src = block[a_ls, b_ls]
            if p == q:
                dst = out[a_ps, b_ps]
                np.add(src, src.T, out=dst)
                np.fill_diagonal(dst, np.diagonal(src))
            else:
                out[a_ps, b_ps] = src
                out[b_ps, a_ps] = src.T
            return
        sub = block[np.ix_(a_locals, b_locals)]
        if p == q:
            res = sub + sub.T
            np.fill_diagonal(res, np.diagonal(sub))
            out[np.ix_(a_positions, b_positions)] = res
        else:
            out[np.ix_(a_positions, b_positions)] = sub
            out[np.ix_(b_positions, a_positions)] = sub.T

    if len(pairs) <= 1 or max_workers <= 1:
        for p, q in pairs:
            blk = reader.read_block(p, q)
            if blk is not None:
                _scatter(p, q, blk)
            del blk
        return out

    workers = min(len(pairs), max_workers)
    max_inflight = max(2, workers * 2)
    pair_iter = iter(pairs)
    inflight: dict[Any, tuple[int, int]] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in range(min(max_inflight, len(pairs))):
            p, q = next(pair_iter)
            inflight[pool.submit(reader.read_block, p, q)] = (p, q)
        while inflight:
            done, _ = wait(inflight, return_when=FIRST_COMPLETED)
            for fut in done:
                p, q = inflight.pop(fut)
                blk = fut.result()
                if blk is not None:
                    _scatter(p, q, blk)
                del blk
                nxt = next(pair_iter, None)
                if nxt is not None:
                    inflight[pool.submit(reader.read_block, nxt[0], nxt[1])] = nxt
    return out


def _apply_flips_inplace(matrix: Any, flip_positions: list[int]) -> Any:
    """Negate rows/cols for swapped-allele variants in-place in 512-row tiles (0 extra O(M^2) allocations)."""
    if not flip_positions:
        return matrix
    import numpy as np

    n = matrix.shape[0]
    signs = np.ones(n, dtype=np.float64)
    signs[np.asarray(flip_positions, dtype=np.intp)] = -1.0
    step = 512
    for i in range(0, n, step):
        ie = min(n, i + step)
        tile = matrix[i:ie, :]
        tile *= signs[i:ie, None]
        tile *= signs[None, :]
    return matrix


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    no_cli = (argv is None and len(sys.argv) == 1)
    default_bm = (
        str(DEFAULT_LOCAL_TEST_BM)
        if (no_cli and DEFAULT_LOCAL_TEST_BM.is_dir())
        else None
    )
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen-json", default=str(DEFAULT_SCREEN_JSON))
    p.add_argument("--locus-index", type=int, default=0, help="0-based index in screen_summary.selected_regions")
    p.add_argument("--region-id", default=None, help="Select locus by region_id instead of --locus-index")
    p.add_argument("--outdir", default=None, help="Output locus folder (defaults to <test_nppb>/loci/<region_id>)")
    p.add_argument("--chain", default=str(DEFAULT_REF_ROOT / "genome/hg38ToHg19.over.chain.gz"))
    p.add_argument("--grch37-fasta", default=str(DEFAULT_REF_ROOT / "genome/GRCh37.fa"))
    p.add_argument("--panukb-dir", default=None)
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR))
    p.add_argument("--bm-path", default=default_bm, help="Optional BlockMatrix path override (defaults to live S3)")
    p.add_argument(
        "--bm-disk-cache",
        default=os.environ.get("GWASFM_BM_CACHE_DIR", ""),
        help="Optional persistent disk cache root for remote BlockMatrix part-* files",
    )
    p.add_argument("--local-ld-panel", default=None, help="Optional offline .npz fixture LD panel")
    p.add_argument("--fetch-workers", type=int, default=4)
    p.add_argument("--block-cache", type=int, default=4)
    p.add_argument("--keep-intermediates", action="store_true")
    return p.parse_args(argv)


def prep_locus(args: argparse.Namespace) -> dict[str, Any]:
    t0 = time.monotonic()
    u0, s0 = current_cpu_times()

    screen_path = Path(args.screen_json).expanduser().resolve()
    if not screen_path.is_file():
        raise FileNotFoundError(
            f"screen_summary.json not found at {screen_path} (run scripts/screen_gwas.py first)."
        )
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    regions = screen.get("selected_regions") or []
    if not regions:
        raise ValueError(f"No active regions in {screen_path}")

    if args.region_id:
        matches_reg = [r for r in regions if r["region_id"] == args.region_id]
        if not matches_reg:
            raise ValueError(f"Region {args.region_id} not found in {screen_path}")
        locus = matches_reg[0]
    else:
        if args.locus_index < 0 or args.locus_index >= len(regions):
            raise IndexError(f"locus_index {args.locus_index} out of range (0..{len(regions)-1})")
        locus = regions[args.locus_index]

    region_id = locus["region_id"]
    region_type = locus["region_type"]
    chrom = normalize_chromosome(locus["CHROM"])
    start = int(locus["start"])
    end = int(locus["end"])

    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else (DEFAULT_LOCUS_ROOT / region_id).resolve()
    )
    outdir.mkdir(parents=True, exist_ok=True)

    vcf = Path(screen["vcf"]).resolve()
    vcf_tags = screen.get("vcf_tags") or inspect_vcf_tags(vcf)
    population = screen["population"]
    gene = screen.get("gene")
    gene_window_bp = int(screen.get("gene_window_bp", 1_000_000))
    chain = Path(args.chain).expanduser().resolve()
    fasta = Path(args.grch37_fasta).expanduser().resolve()
    cache_dir = Path(args.cache_dir).expanduser().resolve()
    panukb_dir = (
        Path(args.panukb_dir).expanduser().resolve()
        if args.panukb_dir
        else default_panukb_dir()
    )
    local_ld_panel = Path(args.local_ld_panel).expanduser().resolve() if args.local_ld_panel else None

    # --- 2. Regional VCF extraction -------------------------------------------
    t_vcf0 = time.monotonic()
    records = query_vcf(
        vcf,
        region=format_vcf_region(vcf_tags, chrom, start, end),
        vcf_tags=vcf_tags,
    )
    vcf_query_seconds = round(time.monotonic() - t_vcf0, 3)

    # --- 3. GRCh38 -> GRCh37 liftover + FASTA-verified normalisation ----------
    t_lift0 = time.monotonic()
    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    with tempfile.TemporaryDirectory(prefix=f"gwasfm-lift-{region_id[:20]}-") as lift_tmp:
        lifted_vcf_path = outdir / "lifted.grch37.vcf" if args.keep_intermediates else Path(lift_tmp) / "lifted.vcf"
        lifted_by_id, liftover_excluded = write_lifted_vcf(records, lifter, lifted_vcf_path)
        norm_tsv_out = outdir / "normalized.grch37.tsv" if args.keep_intermediates else None
        normalized = normalize_lifted_vcf(lifted_vcf_path, fasta, norm_tsv_out) if lifted_by_id else []
    liftover_seconds = round(time.monotonic() - t_lift0, 3)

    for row in normalized:
        row.update({
            key: value
            for key, value in lifted_by_id[row["source_variant_key"]].items()
            if key not in ("LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "ES", "SE", "LP", "SS")
        })

    cross_chrom_excluded = [r for r in normalized if r["LD_CHROM"] != chrom]
    for r in cross_chrom_excluded:
        r["status"] = "excluded"
        r["reason"] = "lifted_to_different_chromosome"
    normalized = [r for r in normalized if r["LD_CHROM"] == chrom]

    if not normalized:
        qc_skip = {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "status": "skipped_no_lifted_variants",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": len(records),
            "variants_lifted": 0,
            "variants_ld_matched": 0,
            "exclusion_reasons": dict(Counter(r["reason"] for r in liftover_excluded + cross_chrom_excluded)),
        }
        write_json(outdir / "locus_prep_qc.json", qc_skip)
        return qc_skip

    target_chrom = chrom
    target_start = max(1, min(r["LD_POS"] for r in normalized) - 100)
    target_end = max(r["LD_POS"] for r in normalized) + 100
    grch37_span = target_end - target_start
    if grch37_span > 9_990_000:
        raise ValueError(
            f"Region {region_id} spans {grch37_span:,} bp on GRCh37, exceeding Pan-UKB's 10 Mb LD band."
        )

    # --- 4. Regional Pan-UKB VariantIndex slice & vectorized allele match -----
    t_idx0 = time.monotonic()
    index_path = ensure_regional_index(
        population, target_chrom, target_start, target_end, cache_dir, SCRIPT_DIR, panukb_dir
    )
    index_cut_seconds = round(time.monotonic() - t_idx0, 3)

    t_m0 = time.monotonic()
    matches = match_panel_variants(index_path, normalized)
    matched, panel_excluded = resolve_panel_matches(normalized, matches)
    match_seconds = round(time.monotonic() - t_m0, 3)

    all_excluded = liftover_excluded + cross_chrom_excluded + panel_excluded
    reason_counts = dict(Counter(r["reason"] for r in all_excluded))

    if not matched:
        qc_skip = {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "target_region_grch37": f"{target_chrom}:{target_start}-{target_end}",
            "status": "skipped_no_panel_matched_variants",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": len(records),
            "variants_lifted": len(normalized),
            "variants_ld_matched": 0,
            "exclusion_reasons": reason_counts,
        }
        write_json(outdir / "locus_prep_qc.json", qc_skip)
        return qc_skip

    matched.sort(key=lambda r: (chromosome_sort_key(r["LD_CHROM"]), r["LD_POS"], r["LD_REF"], r["LD_ALT"]))
    original_by_key = {(r["CHROM"], r["POS"], r["REF"], r["ALT"]): r for r in records}

    matched_rows: list[dict[str, Any]] = []
    for row in matched:
        source = original_by_key.get((row["CHROM"], row["POS"], row["REF"], row["ALT"]))
        if source is None:
            raise ValueError("Could not recover source GWAS variant row.")
        cis_trans_label = (
            "cis"
            if gene is not None
            and source["CHROM"] == gene["chrom"]
            and (gene["start"] - gene_window_bp) <= source["POS"] <= (gene["end"] + gene_window_bp)
            else ("trans" if gene is not None else "na")
        )
        gwas_variant_id = f"{source['CHROM']}:{source['POS']}:{source['REF']}:{source['ALT']}"
        matched_rows.append({
            **source,
            "GWAS_VARIANT_ID": gwas_variant_id,
            "GWAS_CHROM": source["CHROM"],
            "GWAS_POS": source["POS"],
            "GWAS_REF": source["REF"],
            "GWAS_ALT": source["ALT"],
            "CHROM": row["LD_CHROM"],
            "POS": row["LD_POS"],
            "REF": row["LD_REF"],
            "ALT": row["LD_ALT"],
            "LD_CHROM": row["LD_CHROM"],
            "LD_POS": row["LD_POS"],
            "LD_REF": row["LD_REF"],
            "LD_ALT": row["LD_ALT"],
            "LD_VARIANT_ID": row["LD_VARIANT_ID"],
            "PANEL_IDX": row["PANEL_IDX"],
            "PANEL_ALLELE_SWAPPED": row["PANEL_ALLELE_SWAPPED"],
            "EZ": source["ES"] / source["SE"],
            "cis_trans": cis_trans_label,
            "mhc_region": is_mhc_grch38(source["CHROM"], source["POS"]),
            "mhc_core_grc": (
                is_mhc_core_grc_grch38(source["CHROM"], source["POS"])
                or (
                    row["LD_CHROM"] == MHC_CHROM
                    and MHC_CORE_GRC_GRCH37[0] <= row["LD_POS"] <= MHC_CORE_GRC_GRCH37[1]
                )
            ),
        })

    # --- 5. BlockMatrix LD extraction + in-place Fortran .npy write -----------
    t_ld0 = time.monotonic()
    matrix_npy_path = outdir / "matrix.npy"
    lead_idx = max(range(len(matched_rows)), key=lambda idx: (matched_rows[idx]["LP"], -matched_rows[idx]["GWAS_POS"]))

    bm_cache_hits = 0
    bm_cache_misses = 0
    if local_ld_panel is not None:
        import numpy as np

        ld_prefix = outdir / f"panukb_{population.lower()}"
        extract_ld_from_local_panel(local_ld_panel, matched_rows, ld_prefix)
        ld_validation = matrix_to_npy(
            Path(f"{ld_prefix}.npz"),
            Path(f"{ld_prefix}.variants.tsv"),
            matrix_npy_path,
        )
        arr = np.load(matrix_npy_path, mmap_mode="r")
        lead_r2 = (arr[:, lead_idx] ** 2).tolist()
        del arr
        if not args.keep_intermediates:
            Path(f"{ld_prefix}.npz").unlink(missing_ok=True)
            Path(f"{ld_prefix}.variants.tsv").unlink(missing_ok=True)
    else:
        bm_uri = (
            args.bm_path
            if args.bm_path
            else f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.bm"
        )
        global_idxs = [int(r["PANEL_IDX"]) for r in matched_rows]
        reader = open_blockmatrix_reader(
            bm_uri,
            block_cache=min(int(args.block_cache), 2),
            disk_cache_root=args.bm_disk_cache,
        )
        matrix = _assemble_bounded(reader, global_idxs, max_workers=args.fetch_workers)
        bm_cache_hits = int(getattr(reader, "disk_hits", 0))
        bm_cache_misses = int(getattr(reader, "disk_misses", 0))
        del reader
        gc.collect()
        flip_positions = [i for i, r in enumerate(matched_rows) if r["PANEL_ALLELE_SWAPPED"]]
        matrix = _apply_flips_inplace(matrix, flip_positions)

        ids = [
            f"{normalize_chromosome(str(r['LD_CHROM']))}:{int(r['LD_POS'])}:{r['LD_REF']}:{r['LD_ALT']}"
            for r in matched_rows
        ]
        ld_validation = _normalize_ld_matrix_inplace(matrix, ids)
        lead_r2 = (matrix[:, lead_idx] ** 2).tolist()
        _write_symmetric_fortran_npy(matrix, ids, matrix_npy_path)
        del matrix
        gc.collect()
    ld_extract_seconds = round(time.monotonic() - t_ld0, 3)

    for idx, r2_val in enumerate(lead_r2):
        matched_rows[idx]["LD_R2_TO_LEAD"] = r2_val

    matched_variants_path = outdir / "matched_variants.tsv"
    write_tsv(
        matched_variants_path,
        matched_rows,
        [
            "CHROM", "POS", "ID", "REF", "ALT",
            "GWAS_VARIANT_ID", "GWAS_CHROM", "GWAS_POS", "GWAS_REF", "GWAS_ALT",
            "LP", "ES", "SE", "SS", "AF", "AFS", "HM", "EZ", "LD_R2_TO_LEAD",
            "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "LD_VARIANT_ID",
            "PANEL_IDX", "PANEL_ALLELE_SWAPPED", "cis_trans",
            "mhc_region", "mhc_core_grc",
        ],
    )

    n_swapped = sum(1 for r in matched_rows if r["PANEL_ALLELE_SWAPPED"])
    lead_matched = matched_rows[lead_idx]
    ss_vals = sorted(float(r["SS"]) for r in matched_rows if r.get("SS") is not None)
    ss_min = ss_vals[0] if ss_vals else None
    ss_max = ss_vals[-1] if ss_vals else None
    ss_median = ss_vals[len(ss_vals) // 2] if ss_vals else None
    ss_ratio = round(ss_min / ss_max, 4) if (ss_min is not None and ss_max and ss_max > 0) else None
    ss_het_warn = bool(ss_ratio is not None and ss_ratio < 0.80)

    u1, s1 = current_cpu_times()
    wall_seconds = round(time.monotonic() - t0, 3)
    npy_mb = round(matrix_npy_path.stat().st_size / (1024.0 * 1024.0), 2)

    prep_qc: dict[str, Any] = {
        "region_id": region_id,
        "region_type": region_type,
        "status": "prepared",
        "source_region_grch38": f"{chrom}:{start}-{end}",
        "target_region_grch37": f"{target_chrom}:{target_start}-{target_end}",
        "width_bp_grch38": end - start,
        "width_bp_grch37": grch37_span,
        "mhc_region": locus["mhc_region"],
        "mhc_core_grc": any(r["mhc_core_grc"] for r in matched_rows),
        "cis_overlap": locus["cis_overlap"],
        "split_at_panukb_radius": locus.get("split_at_panukb_radius", False),
        "extended_by_adjacent_lead": locus.get("extended_by_adjacent_lead", False),
        "significant_variants_in_region": locus["significant_variant_count"],
        "variants_in_source_region": len(records),
        "variants_lifted": len(normalized),
        "variants_ld_matched": len(matched_rows),
        "variants_allele_swapped": n_swapped,
        "variants_excluded": len(records) - len(matched_rows),
        "exclusion_reasons": reason_counts,
        "ss_min": ss_min,
        "ss_median": ss_median,
        "ss_max": ss_max,
        "ss_min_max_ratio": ss_ratio,
        "sample_size_heterogeneity_warning": ss_het_warn,
        "lead_variant_id": lead_matched["ID"],
        "lead_gwas_variant_id": lead_matched["GWAS_VARIANT_ID"],
        "lead_pos_grch38": lead_matched["GWAS_POS"],
        "lead_lp": lead_matched["LP"],
        "index_cut_seconds": index_cut_seconds,
        "wall_seconds": wall_seconds,
        "ld_validation": ld_validation,
        "bm_source": "local_override" if args.bm_path else ("local_npz" if local_ld_panel else "live_s3"),
        "bm_cache_hits": bm_cache_hits,
        "bm_cache_misses": bm_cache_misses,
        "scratch_matrix_mb": npy_mb,
        "timing": {
            "vcf_query_seconds": vcf_query_seconds,
            "liftover_norm_seconds": liftover_seconds,
            "index_cut_seconds": index_cut_seconds,
            "variant_match_seconds": match_seconds,
            "ld_extract_normalize_seconds": ld_extract_seconds,
            "wall_seconds": wall_seconds,
            "user_seconds": round(u1 - u0, 3),
            "sys_seconds": round(s1 - s0, 3),
            "peak_rss_mb": current_peak_rss_mb(),
        },
    }
    write_json(outdir / "locus_prep_qc.json", prep_qc)
    print(
        f"[prep_locus] {region_id}: {len(matched_rows)}/{len(records)} matched | "
        f"matrix={npy_mb:.1f} MB | bm_cache={bm_cache_hits}h/{bm_cache_misses}m | "
        f"ld_extract={ld_extract_seconds:.2f}s, total={wall_seconds:.2f}s, "
        f"peak RSS {prep_qc['timing']['peak_rss_mb']:.1f} MB -> {outdir}"
    )
    return prep_qc


if __name__ == "__main__":
    prep_locus(parse_args())
