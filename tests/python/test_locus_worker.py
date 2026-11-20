"""Unit tests for coordinated node-local /tmp space reservation and spillover."""

from __future__ import annotations

import os
import sys
import time
import unittest
from pathlib import Path
from unittest import mock

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import locus_worker  # noqa: E402
import screen_gwas  # noqa: E402


class LocusWorkerScratchTests(unittest.TestCase):
    def test_automatic_locus_memory_scaling_by_variants_and_width(self) -> None:
        # Sparse HapMap 3 Mb locus (e.g. Locke BMI ~3.4k variants) -> 16 GB
        self.assertEqual(screen_gwas.estimate_locus_memory_gb(3_400, 3_000_000, "lead_locus"), 16)
        # Standard UKB-PPP 3 Mb locus (~20k variants) -> 16 GB
        self.assertEqual(screen_gwas.estimate_locus_memory_gb(20_000, 3_000_000, "lead_locus"), 16)
        # Dense UKB-PPP 3 Mb locus (~24k variants) -> 24 GB
        self.assertEqual(screen_gwas.estimate_locus_memory_gb(24_000, 3_000_000, "lead_locus"), 24)
        # Wide 9.33 Mb merged UKB-PPP locus (~65k variants, ~31.5 GiB matrix) -> 128 GB
        self.assertEqual(screen_gwas.estimate_locus_memory_gb(65_000, 9_333_934, "lead_locus"), 128)
        # Extended 8.50 Mb UKB-PPP MHC locus (~80k variants, ~47.7 GiB matrix) -> 160 GB
        self.assertEqual(screen_gwas.estimate_locus_memory_gb(80_000, 8_500_000, "mhc"), 160)
        # Sparse HapMap MHC locus (~10k variants) stays small -> 16 GB
        self.assertEqual(screen_gwas.estimate_locus_memory_gb(10_000, 8_500_000, "mhc"), 16)
    def test_concurrent_reservations_spill_cleanly_and_reclaim_on_release(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            node_tmp = base / "node_tmp"
            spill_root = base / "spill_scratch"
            node_tmp.mkdir()
            spill_root.mkdir()

            # Simulate a 20 GiB free /tmp with 5.5 GiB per-worker reserve + 4 GiB headroom -> max 2 workers on /tmp
            fake_usage = os.statvfs_result if hasattr(os, "statvfs_result") else None
            _ = fake_usage
            gib = 1024 ** 3
            with mock.patch("locus_worker.shutil.disk_usage") as m_du:
                m_du.return_value = mock.Mock(total=24 * gib, used=4 * gib, free=20 * gib)

                c1 = locus_worker.allocate_locus_workdir(
                    "chr1_1_100",
                    spill_root,
                    reserve_bytes=int(5.5 * gib),
                    min_headroom_bytes=4 * gib,
                    tmp_candidates=[node_tmp],
                )
                c2 = locus_worker.allocate_locus_workdir(
                    "chr1_2_200",
                    spill_root,
                    reserve_bytes=int(5.5 * gib),
                    min_headroom_bytes=4 * gib,
                    tmp_candidates=[node_tmp],
                )
                # 3rd concurrent claim: 20 - 11 - 5.5 = 3.5 GiB < 4 GiB headroom -> spills to spill_root!
                c3 = locus_worker.allocate_locus_workdir(
                    "chr1_3_300",
                    spill_root,
                    reserve_bytes=int(5.5 * gib),
                    min_headroom_bytes=4 * gib,
                    tmp_candidates=[node_tmp],
                )

                self.assertTrue(c1.using_node_tmp)
                self.assertEqual(c1.work_dir.parent, node_tmp)
                self.assertTrue(c2.using_node_tmp)
                self.assertEqual(c2.work_dir.parent, node_tmp)
                self.assertFalse(c3.using_node_tmp)
                self.assertEqual(c3.work_dir.parent, spill_root)

                # Releasing c1 frees its 5.5 GiB slot on node_tmp immediately
                locus_worker.release_locus_workdir(c1)
                self.assertFalse(c1.work_dir.exists())

                c4 = locus_worker.allocate_locus_workdir(
                    "chr1_4_400",
                    spill_root,
                    reserve_bytes=int(5.5 * gib),
                    min_headroom_bytes=4 * gib,
                    tmp_candidates=[node_tmp],
                )
                self.assertTrue(c4.using_node_tmp)
                self.assertEqual(c4.work_dir.parent, node_tmp)

                locus_worker.release_locus_workdir(c2)
                locus_worker.release_locus_workdir(c3)
                locus_worker.release_locus_workdir(c4)

    def test_sweeps_dead_pid_and_legacy_orphan_directories_while_preserving_other_hosts(self) -> None:
        import socket
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            node_tmp = Path(td)
            my_host = socket.gethostname()

            dead_dir = node_tmp / "gwasfm_chr2_dead_abc"
            dead_dir.mkdir()
            (dead_dir / ".owner_pid").write_text(f"{my_host}:99999999\n", encoding="utf-8")

            other_host_dir = node_tmp / "gwasfm_chr2_remote_node"
            other_host_dir.mkdir()
            (other_host_dir / ".owner_pid").write_text("cn99_other_host:99999999\n", encoding="utf-8")

            legacy_dir = node_tmp / "gwasfm_chr3_legacy_xyz"
            legacy_dir.mkdir()
            old_ts = time.time() - 7200.0
            os.utime(legacy_dir, (old_ts, old_ts))

            live_dir = node_tmp / "gwasfm_chr4_live_def"
            live_dir.mkdir()
            (live_dir / ".owner_pid").write_text(f"{my_host}:{os.getpid()}\n", encoding="utf-8")

            removed = locus_worker._sweep_stale_tmp_dirs(node_tmp, os.getuid())
            self.assertEqual(removed, 2)
            self.assertFalse(dead_dir.exists())
            self.assertFalse(legacy_dir.exists())
            self.assertTrue(other_host_dir.exists())
            self.assertTrue(live_dir.exists())

    def test_home_and_shared_tmpdir_excluded_from_node_tmp_candidates(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            fake_home = Path(td) / "home_leem"
            home_tmp = fake_home / "tmp"
            shared_spill = Path(td) / "shared_scratch_tmp"
            slurm_local = Path(td) / "slurm_local_tmp"
            home_tmp.mkdir(parents=True)
            shared_spill.mkdir(parents=True)
            slurm_local.mkdir(parents=True)

            roots = locus_worker.candidate_node_tmp_roots(
                environ={"SLURM_TMPDIR": str(slurm_local), "TMPDIR": str(shared_spill)},
                home_dir=fake_home,
                exclude_roots=[shared_spill],
            )
            resolved_strs = [str(p) for p in roots]
            self.assertNotIn(str(home_tmp.resolve()), resolved_strs)
            self.assertNotIn(str(shared_spill.resolve()), resolved_strs)
            self.assertIn(str(slurm_local.resolve()), resolved_strs)

    def test_enospc_detection(self) -> None:
        self.assertTrue(
            locus_worker.is_enospc_failure(
                "OSError: Not enough free space to write 3260765768 bytes after offset 128"
            )
        )
        self.assertTrue(locus_worker.is_enospc_failure("OSError: [Errno 28] No space left on device"))
        self.assertFalse(locus_worker.is_enospc_failure("ValueError: invalid allele"))

    def test_explicit_temp_workdir_stages_inside_requested_dir_and_cleans_up(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            explicit_root = Path(td) / "scratch_root" / "tmp"
            claim = locus_worker.allocate_explicit_temp_workdir("chr6_1_2", explicit_root)
            self.assertFalse(claim.using_node_tmp)
            self.assertTrue(claim.cleanup_on_exit)
            self.assertEqual(claim.work_dir.parent, explicit_root)
            self.assertTrue(claim.work_dir.is_dir())
            self.assertTrue((claim.work_dir / ".owner_pid").is_file())

            locus_worker.release_locus_workdir(claim)
            self.assertFalse(claim.work_dir.exists())

    def test_final_attempt_failure_records_failed_locus_qc_for_isolation(self) -> None:
        import json
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "loci"
            loc_dir = root / "chr6_25500000_34000000_mhc"
            out_file = Path(td) / ".scratch_loci_root"
            rc = locus_worker._record_final_locus_failure(
                loc_dir,
                "chr6_25500000_34000000_mhc",
                "mhc",
                "prep_locus failed: simulated OOM",
                str(out_file),
                root,
            )
            self.assertEqual(rc, 0)
            qc = json.loads((loc_dir / "locus_prep_qc.json").read_text(encoding="utf-8"))
            self.assertEqual(qc["status"], "failed_locus")
            self.assertIn("simulated OOM", qc["error"])
            self.assertEqual(out_file.read_text(encoding="utf-8"), str(root))

    def test_disk_cached_blockmatrix_reader_caches_metadata_and_parts_atomically(self) -> None:
        import io
        import json
        import tempfile
        import prep_locus
        import numpy as np

        raw_block = np.eye(4, dtype=np.float64).tobytes()
        meta_bytes = json.dumps({
            "blockSize": 4,
            "nRows": 4,
            "nCols": 4,
            "maybeFiltered": [0],
            "partFiles": ["part-0-test"],
        }).encode("utf-8")

        fetch_counts = {"meta": 0, "part": 0}

        class FakeFS:
            def open(self, path: str, mode: str = "rb"):
                if path.endswith("/metadata.json"):
                    fetch_counts["meta"] += 1
                    return io.BytesIO(meta_bytes)
                if path.endswith("/parts/part-0-test"):
                    fetch_counts["part"] += 1
                    return io.BytesIO(raw_block)
                raise FileNotFoundError(path)

        with tempfile.TemporaryDirectory() as td:
            cache_dir = Path(td) / "bm_cache"
            with (
                mock.patch("ldcov.io.fs_utils.resolve_filesystem", return_value=(FakeFS(), "")),
                mock.patch(
                    "ldcov.io.blockmatrix.block_codec.decode_block",
                    side_effect=lambda b: np.frombuffer(b, dtype=np.float64).reshape(4, 4).copy(),
                ),
            ):
                r1 = prep_locus.open_blockmatrix_reader(
                    "s3://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm",
                    block_cache=2,
                    disk_cache_root=str(cache_dir),
                )
                b1 = r1.read_block(0, 0)
                self.assertEqual(r1.disk_misses, 1)
                self.assertEqual(r1.disk_hits, 0)
                self.assertEqual(fetch_counts, {"meta": 1, "part": 1})
                np.testing.assert_allclose(b1, np.eye(4))

                # Second reader hits local disk cache with zero remote FS calls!
                r2 = prep_locus.open_blockmatrix_reader(
                    "s3://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm",
                    block_cache=2,
                    disk_cache_root=str(cache_dir),
                )
                b2 = r2.read_block(0, 0)
                self.assertEqual(r2.disk_hits, 1)
                self.assertEqual(r2.disk_misses, 0)
                self.assertEqual(fetch_counts, {"meta": 1, "part": 1})
                np.testing.assert_allclose(b2, np.eye(4))

    def test_per_locus_scratch_keying_completeness_skip_and_migration(self) -> None:
        import json
        import tempfile

        methods = ["p_ld", "susie_rss", "susie_rrss", "carma", "polyfun"]
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            scratch_root = base / "scratch"
            outdir = base / "finemapping"
            rows_dir = base / "rows"
            old_run_scratch = scratch_root / "1005978"
            rows_dir.mkdir(parents=True)

            row = {
                "study_id": "sun_2023_37794186",
                "gwas_id": "P00749_OID21124",
                "settings_hash": "a1b2c3d4e5f6",
                "outdir": str(outdir),
                "methods": methods,
                "force": False,
                "keep_intermediates": False,
            }
            row_json = rows_dir / "row_00027_sun_2023_37794186__P00749_OID21124.json"
            row_json.write_text(json.dumps(row), encoding="utf-8")

            screen_json = base / "screen_summary.json"
            screen_json.write_text(
                json.dumps({
                    "reference_sample_size": 420542,
                    "selected_regions": [{"region_id": "chr9_1_100", "vcf_variants_in_region": 22000}],
                }),
                encoding="utf-8",
            )

            # Populate 1 complete locus + 1 incomplete locus in legacy run scratch
            old_loci = old_run_scratch / row_json.stem / "loci"
            good_loc = old_loci / "chr9_1_100"
            good_loc.mkdir(parents=True)
            (good_loc / "locus_prep_qc.json").write_text(json.dumps({"status": "prepared"}), encoding="utf-8")
            for m in methods:
                md = good_loc / "methods" / m
                md.mkdir(parents=True)
                (md / "region_result.rds").write_bytes(b"RDSDATA")
                (md / "method_qc.json").write_text(json.dumps({"method_id": m}), encoding="utf-8")

            bad_loc = old_loci / "chr6_25500000_34000000_mhc"
            bad_loc.mkdir(parents=True)
            (bad_loc / "locus_prep_qc.json").write_text(json.dumps({"status": "prepared"}), encoding="utf-8")
            # Only 2 of 5 methods present -> incomplete!
            for m in ("p_ld", "carma"):
                md = bad_loc / "methods" / m
                md.mkdir(parents=True)
                (md / "region_result.rds").write_bytes(b"RDSDATA")
                (md / "method_qc.json").write_text(json.dumps({"method_id": m}), encoding="utf-8")

            self.assertTrue(locus_worker.locus_outputs_are_complete(good_loc, methods))
            self.assertFalse(locus_worker.locus_outputs_are_complete(bad_loc, methods))

            stats = locus_worker.migrate_run_scratch(old_run_scratch, rows_dir, scratch_root)
            self.assertEqual(stats, {"migrated_loci": 1, "gwas_with_migrated_loci": 1})

            dst_loci = locus_worker.gwas_scratch_loci_dir(scratch_root, row)
            self.assertTrue(locus_worker.locus_outputs_are_complete(dst_loci / "chr9_1_100", methods))
            self.assertFalse((dst_loci / "chr6_25500000_34000000_mhc").exists())

            # Invoking run_locus_worker on the migrated complete locus skips subprocess execution completely
            out_root_file = base / ".scratch_loci_root"
            args = locus_worker.parse_args([
                "--project-dir", str(SCRIPTS_DIR.parent),
                "--row-json", str(row_json),
                "--screen-json", str(screen_json),
                "--locus-index", "0",
                "--region-id", "chr9_1_100",
                "--region-type", "lead_locus",
                "--scratch-root", str(scratch_root),
                "--out-root-file", str(out_root_file),
            ])
            with mock.patch("locus_worker._run_and_log") as m_run:
                rc = locus_worker.run_locus_worker(args)
                self.assertEqual(rc, 0)
                m_run.assert_not_called()
            self.assertEqual(out_root_file.read_text(encoding="utf-8"), str(dst_loci))

    def test_dynamic_locus_reserve_bytes_and_failed_gwas_preserves_good_loci(self) -> None:
        import json
        import tempfile
        import finalise_gwas

        screen = {
            "selected_regions": [
                {"region_id": "chr1_1_10", "vcf_variants_in_region": 3400},
                {"region_id": "chr9_1_10", "vcf_variants_in_region": 22000},
                {"region_id": "chr6_25500000_34000000_mhc", "vcf_variants_in_region": 80000},
            ]
        }
        r_sparse = locus_worker.estimate_locus_reserve_bytes(screen, 0, "lead_locus")
        r_std = locus_worker.estimate_locus_reserve_bytes(screen, 1, "lead_locus")
        r_mhc = locus_worker.estimate_locus_reserve_bytes(screen, 2, "mhc")
        self.assertEqual(r_sparse, int(1.5 * (1024 ** 3)))
        self.assertGreater(r_std, int(4.0 * (1024 ** 3)))
        self.assertLess(r_std, int(6.0 * (1024 ** 3)))
        self.assertGreater(r_mhc, int(50.0 * (1024 ** 3)))

        # Verify finalise_gwas removes ONLY failed locus directory while preserving completed sibling loci
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            loci_root = base / "scratch" / "by_gwas" / "s__g__h" / "loci"
            good_loc = loci_root / "chr1_1_10"
            (good_loc / "methods" / "p_ld").mkdir(parents=True)
            (good_loc / "locus_prep_qc.json").write_text(
                json.dumps({"region_id": "chr1_1_10", "region_type": "lead_locus", "status": "prepared"}),
                encoding="utf-8",
            )
            (good_loc / "methods" / "p_ld" / "region_result.rds").write_bytes(b"RDS")
            (good_loc / "methods" / "p_ld" / "method_qc.json").write_text(
                json.dumps({"method_id": "p_ld"}), encoding="utf-8"
            )
            bad_loc = loci_root / "chr6_25500000_34000000_mhc"
            bad_loc.mkdir(parents=True)
            (bad_loc / "locus_prep_qc.json").write_text(
                json.dumps({
                    "region_id": "chr6_25500000_34000000_mhc",
                    "region_type": "mhc",
                    "status": "failed_locus",
                    "error": "OOM",
                }),
                encoding="utf-8",
            )
            screen_file = base / "screen_summary.json"
            screen_file.write_text(
                json.dumps({
                    "study_id": "s",
                    "gwas_id": "g",
                    "population": "EUR",
                    "ancestry_ld_status": "matched",
                    "reference_sample_size": 420542,
                    "significance_lp": 7.3,
                    "region_flank_bp": 1500000,
                    "genomewide_significant_variants": 10,
                    "candidate_regions_total": 2,
                    "selected_regions_count": 2,
                    "selected_regions": [
                        {"region_id": "chr1_1_10", "region_type": "lead_locus"},
                        {"region_id": "chr6_25500000_34000000_mhc", "region_type": "mhc"},
                    ],
                }),
                encoding="utf-8",
            )
            f_args = finalise_gwas.parse_args([
                "--screen-json", str(screen_file),
                "--loci-dir", str(loci_root),
                "--output-root", str(base / "out"),
                "--methods", "p_ld",
                "--cleanup-loci-scratch",
                "--force",
            ])
            with self.assertRaises(RuntimeError):
                finalise_gwas.finalise_gwas(f_args)
            self.assertTrue(good_loc.is_dir())
            self.assertFalse(bad_loc.exists())

    def test_partial_locus_resumes_only_missing_methods(self) -> None:
        import json
        import tempfile

        methods = ["p_ld", "susie_rss", "susie_rrss", "carma", "polyfun"]
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            scratch_root = base / "scratch"
            outdir = base / "finemapping"
            row = {
                "study_id": "sun_2023_37794186",
                "gwas_id": "O15467_OID20334",
                "settings_hash": "ba2586a3922cdc95",
                "outdir": str(outdir),
                "methods": methods,
                "force": False,
                "keep_intermediates": False,
                "chain": "/tmp/chain",
                "grch37_fasta": "/tmp/fa",
                "cache_dir": "/tmp/cache",
            }
            row_json = base / "row.json"
            row_json.write_text(json.dumps(row), encoding="utf-8")
            screen_json = base / "screen.json"
            screen_json.write_text(
                json.dumps({
                    "reference_sample_size": 420542,
                    "selected_regions": [{"region_id": "chr6_61385886_67594203", "vcf_variants_in_region": 51085}],
                }),
                encoding="utf-8",
            )

            loc_dir = locus_worker.gwas_scratch_loci_dir(scratch_root, row) / "chr6_61385886_67594203"
            loc_dir.mkdir(parents=True)
            (loc_dir / "locus_prep_qc.json").write_text(
                json.dumps({"status": "prepared", "variants_ld_matched": 51085}),
                encoding="utf-8",
            )
            # Pre-populate 4 of 5 methods complete (only polyfun missing!)
            for m in ("p_ld", "carma", "susie_rss", "susie_rrss"):
                md = loc_dir / "methods" / m
                md.mkdir(parents=True)
                (md / "region_result.rds").write_bytes(b"SAVED_RDS")
                (md / "method_qc.json").write_text(json.dumps({"method_id": m}), encoding="utf-8")

            recorded_cmds: list[list[str]] = []

            def fake_run(cmd: list[str], step_stem: str, gwas_log_dir: Path, sub_env: dict[str, str]) -> None:
                recorded_cmds.append(cmd)
                if "prep_locus.py" in " ".join(cmd):
                    wd = Path(cmd[cmd.index("--outdir") + 1])
                    (wd / "locus_prep_qc.json").write_text(
                        json.dumps({"status": "prepared", "variants_ld_matched": 51085}),
                        encoding="utf-8",
                    )
                    (wd / "matched_variants.tsv").write_text("CHROM\n", encoding="utf-8")
                    (wd / "matrix.npy").write_bytes(b"NPY")

            args = locus_worker.parse_args([
                "--project-dir", str(SCRIPTS_DIR.parent),
                "--row-json", str(row_json),
                "--screen-json", str(screen_json),
                "--locus-index", "0",
                "--region-id", "chr6_61385886_67594203",
                "--region-type", "lead_locus",
                "--scratch-root", str(scratch_root),
                "--out-root-file", str(base / ".out_root"),
            ])
            with mock.patch("locus_worker._run_and_log", side_effect=fake_run):
                rc = locus_worker.run_locus_worker(args)
            self.assertEqual(rc, 0)
            self.assertEqual(len(recorded_cmds), 2)
            r_cmd = recorded_cmds[1]
            self.assertEqual(r_cmd[r_cmd.index("--methods") + 1], "polyfun")
            for m in ("p_ld", "carma", "susie_rss", "susie_rrss"):
                self.assertTrue((loc_dir / "methods" / m / "region_result.rds").is_file())

    def test_case_control_vcf_header_sample_size_without_format_ss(self) -> None:
        import run_finemapping

        sample_hdr = (
            '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect">\n'
            '##FORMAT=<ID=SE,Number=A,Type=Float,Description="SE">\n'
            '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10P">\n'
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">\n'
            '##SAMPLE=<ID=GECCO_colon_all_noUKB,TotalVariants=14131476,TotalControls=38626,'
            'TotalCases=27817,StudyType=CaseControl,Phenotype="colon cancer, both sexes",Ancestry=ALL>\n'
        )
        parsed = run_finemapping._parse_vcf_sample_header(sample_hdr)
        self.assertEqual(parsed["SAMPLE_STUDY_TYPE"], "CaseControl")
        self.assertEqual(parsed["SAMPLE_TOTAL_CASES"], 27817)
        self.assertEqual(parsed["SAMPLE_TOTAL_CONTROLS"], 38626)
        expected_neff = round(4.0 / (1.0 / 27817.0 + 1.0 / 38626.0), 2)
        self.assertEqual(parsed["HEADER_EFFECTIVE_N"], expected_neff)
        self.assertEqual(parsed["HEADER_EFFECTIVE_N_SOURCE"], "vcf_sample_header_neff_4_ncases_ncontrols")

        import shutil
        import subprocess
        import tempfile
        if shutil.which("bcftools"):
            with tempfile.TemporaryDirectory() as td:
                vcf_txt = Path(td) / "cc.vcf"
                vcf_gz = Path(td) / "cc.vcf.gz"
                vcf_txt.write_text(
                    "##fileformat=VCFv4.2\n"
                    "##contig=<ID=1,length=248956422>\n"
                    + sample_hdr
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGECCO_colon_all_noUKB\n"
                    + "1\t1000\trs1\tA\tG\t.\t.\t.\tES:SE:LP:AF\t0.12:0.02:8.95:0.35\n",
                    encoding="utf-8",
                )
                subprocess.run(["bcftools", "view", "-Oz", "-o", str(vcf_gz), str(vcf_txt)], check=True)
                subprocess.run(["bcftools", "index", "-t", str(vcf_gz)], check=True)
                tags = run_finemapping.inspect_vcf_tags(vcf_gz)
                self.assertFalse(tags["FORMAT_SS"])
                self.assertEqual(tags["HEADER_EFFECTIVE_N"], expected_neff)
                recs_scan = run_finemapping.query_vcf(vcf_gz, minimum_lp=7.3, vcf_tags=tags)
                recs_reg = run_finemapping.query_vcf(vcf_gz, region="1:500-1500", vcf_tags=tags)
                self.assertEqual(len(recs_scan), 1)
                self.assertEqual(len(recs_reg), 1)
                self.assertEqual(recs_scan[0]["SS"], expected_neff)
                self.assertEqual(recs_reg[0]["SS"], expected_neff)

        # Also verify CaseControl without TotalCases/TotalControls estimates Neff from SE + AF
        hdr_no_n = (
            '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect">\n'
            '##FORMAT=<ID=SE,Number=A,Type=Float,Description="SE">\n'
            '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10P">\n'
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">\n'
            '##SAMPLE=<ID=GECCO_colon_male_noUKB,StudyType=CaseControl,Ancestry=ALL>\n'
        )
        parsed_no_n = run_finemapping._parse_vcf_sample_header(hdr_no_n)
        self.assertNotIn("HEADER_EFFECTIVE_N", parsed_no_n)
        # True Neff = 32000, AF = 0.25 -> SE = sqrt(2 / (0.25 * 0.75 * 32000))
        true_neff = 32000.0
        se_val = (2.0 / (0.25 * 0.75 * true_neff)) ** 0.5
        est = run_finemapping.estimate_case_control_effective_n(
            [{"AF": 0.25, "SE": se_val}],
            parsed_no_n,
        )
        self.assertAlmostEqual(est, true_neff, places=1)
        # Verify /dev/shm is excluded from automatic node temp candidates
        cand_strs = [str(p) for p in locus_worker.candidate_node_tmp_roots()]
        self.assertNotIn("/dev/shm", cand_strs)

    def test_bounded_assemble_and_inplace_flips_match_ldcov_reference(self) -> None:
        import numpy as np
        import prep_locus
        from ldcov.io.blockmatrix.partitioner import GridPartitioner
        from ldcov.ld_bm.extract import _apply_flips, _assemble

        rng = np.random.default_rng(123)
        full = rng.normal(size=(16, 16))
        full = np.corrcoef(full)
        bs = 4
        blocks: dict[tuple[int, int], np.ndarray] = {}
        filtered: list[int] = []
        n_grid = 16 // bs
        for i in range(n_grid):
            for j in range(i, n_grid):
                blk = full[i * bs:(i + 1) * bs, j * bs:(j + 1) * bs].copy()
                if i == j:
                    blk = np.triu(blk)
                blocks[(i, j)] = blk
                filtered.append(i + j * n_grid)

        class FakeBMReader:
            block_size = bs
            partitioner = GridPartitioner(bs, 16, 16, sorted(filtered))

            def read_block(self, i: int, j: int) -> np.ndarray:
                return blocks[(i, j)].copy()

        for idxs in ([0, 1, 2, 3, 4, 5, 8, 9, 12, 15], list(range(16))):
            ref = _assemble(FakeBMReader(), idxs, max_workers=3)
            got = prep_locus._assemble_bounded(FakeBMReader(), idxs, max_workers=3)
            np.testing.assert_allclose(got, ref, rtol=0, atol=1e-15)
            flips = [0, 3, 5, len(idxs) - 1]
            ref_f = _apply_flips(ref.copy(), flips)
            got_f = prep_locus._apply_flips_inplace(got.copy(), flips)
            np.testing.assert_allclose(got_f, ref_f, rtol=0, atol=1e-15)


if __name__ == "__main__":
    unittest.main()
