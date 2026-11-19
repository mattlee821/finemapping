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

    def test_sweeps_dead_pid_and_legacy_orphan_directories(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            node_tmp = Path(td)
            dead_dir = node_tmp / "gwasfm_chr2_dead_abc"
            dead_dir.mkdir()
            (dead_dir / ".owner_pid").write_text("99999999\n", encoding="utf-8")

            legacy_dir = node_tmp / "gwasfm_chr3_legacy_xyz"
            legacy_dir.mkdir()
            old_ts = time.time() - 600.0
            os.utime(legacy_dir, (old_ts, old_ts))

            live_dir = node_tmp / "gwasfm_chr4_live_def"
            live_dir.mkdir()
            (live_dir / ".owner_pid").write_text(f"{os.getpid()}\n", encoding="utf-8")

            removed = locus_worker._sweep_stale_tmp_dirs(node_tmp, os.getuid())
            self.assertEqual(removed, 2)
            self.assertFalse(dead_dir.exists())
            self.assertFalse(legacy_dir.exists())
            self.assertTrue(live_dir.exists())

    def test_home_directory_excluded_from_node_tmp_candidates(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            fake_home = Path(td) / "home_leem"
            home_tmp = fake_home / "tmp"
            real_tmp = Path(td) / "local_tmp"
            home_tmp.mkdir(parents=True)
            real_tmp.mkdir(parents=True)

            roots = locus_worker.candidate_node_tmp_roots(
                environ={"SLURM_TMPDIR": str(home_tmp), "TMPDIR": str(real_tmp)},
                home_dir=fake_home,
            )
            resolved_strs = [str(p) for p in roots]
            self.assertNotIn(str(home_tmp.resolve()), resolved_strs)
            self.assertIn(str(real_tmp.resolve()), resolved_strs)

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


if __name__ == "__main__":
    unittest.main()
