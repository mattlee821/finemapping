import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from run_finemapping import (
    MHC_END_GRCH38,
    MHC_START_GRCH38,
    REFERENCE_SAMPLES,
    cluster_significant,
    cut_regional_index,
    ensure_regional_index,
    greedy_kanai_windows,
    installed_population_index,
    is_mhc_core_grc_grch38,
    is_mhc_grch38,
    locus_window,
    map_variant,
    matrix_to_npy,
    matrix_to_tsv,
    r_runtime_environment,
    resolve_panel_matches,
    resolve_reference_sample_size,
    select_mapped_loci,
    write_lifted_vcf,
)


class StubLifter:
    def __init__(self, mapping):
        self.mapping = mapping

    def convert_coordinate(self, chromosome, position):
        return self.mapping.get((chromosome, position), [])


class FineMappingTests(unittest.TestCase):
    def test_all_six_panukb_reference_sample_sizes_are_wired(self):
        self.assertEqual(
            set(REFERENCE_SAMPLES),
            {"AFR", "AMR", "CSA", "EAS", "EUR", "MID"},
        )
        for pop in REFERENCE_SAMPLES:
            n = resolve_reference_sample_size(pop, None)
            self.assertIsNotNone(n)
            self.assertGreater(n, 500)
        self.assertEqual(REFERENCE_SAMPLES["EUR"], 420_542)
        self.assertEqual(REFERENCE_SAMPLES["AFR"], 6_637)

    def test_greedy_kanai_windows_merges_overlapping_3mb_windows_and_clamps_mhc(self):
        records = [
            # Two leads 2.2 Mb apart -> ±1.5 Mb windows [3.5M, 6.5M] and [5.7M, 8.7M] merge!
            {"CHROM": "1", "POS": 5_000_000, "ID": "rs1", "LP": 15.0},
            {"CHROM": "1", "POS": 5_400_000, "ID": "rs_masked", "LP": 9.0},
            {"CHROM": "1", "POS": 7_200_000, "ID": "rs2", "LP": 11.0},
            # Lead just upstream of MHC (25.5 Mb): clamped at 25,499,999
            {"CHROM": "6", "POS": 24_800_000, "ID": "rs_pre_mhc", "LP": 8.5},
            # Variant inside MHC is skipped by greedy_kanai_windows (handled as single MHC locus)
            {"CHROM": "6", "POS": 31_000_000, "ID": "rs_in_mhc", "LP": 25.0},
        ]
        windows = greedy_kanai_windows(records, flank_bp=1_500_000, max_region_width_bp=9_500_000)
        self.assertEqual(len(windows), 2)
        chr1 = next(w for w in windows if w["CHROM"] == "1")
        self.assertEqual((chr1["start"], chr1["end"]), (3_500_000, 8_700_000))
        self.assertEqual(chr1["lead"]["ID"], "rs1")
        self.assertFalse(chr1["split_at_panukb_radius"])

        chr6 = next(w for w in windows if w["CHROM"] == "6")
        self.assertEqual(chr6["end"], MHC_START_GRCH38 - 1)
        self.assertTrue(is_mhc_grch38("6", 31_000_000))
        self.assertTrue(is_mhc_core_grc_grch38("6", 31_000_000))
        self.assertFalse(is_mhc_core_grc_grch38("6", 26_000_000))
        self.assertTrue(is_mhc_grch38("6", 26_000_000))

    def test_greedy_kanai_splits_hyper_polygenic_chain_above_max_width(self):
        # Chain 5 leads spaced 2.5 Mb apart -> raw merged span = 13.0 Mb (> 9.5 Mb cap)
        records = [
            {"CHROM": "2", "POS": 2_000_000 + i * 2_500_000, "ID": f"rs{i}", "LP": 12.0 - i}
            for i in range(5)
        ]
        windows = greedy_kanai_windows(records, flank_bp=1_500_000, max_region_width_bp=9_500_000)
        self.assertGreaterEqual(len(windows), 2)
        for w in windows:
            self.assertLessEqual(w["end"] - w["start"], 9_500_000)
            self.assertTrue(w["split_at_panukb_radius"])

    def test_clusters_significant_variants_by_chromosome_and_gap(self):
        records = [
            {"CHROM": "2", "POS": 300, "LP": 9},
            {"CHROM": "1", "POS": 101, "LP": 8},
            {"CHROM": "1", "POS": 100, "LP": 10},
            {"CHROM": "1", "POS": 500, "LP": None},
            {"CHROM": "1", "POS": 201, "LP": 7},
        ]

        clusters = cluster_significant(records, merge_distance_bp=100)

        self.assertEqual(
            [[(row["CHROM"], row["POS"]) for row in cluster] for cluster in clusters],
            [[("1", 100), ("1", 101), ("1", 201)], [("2", 300)]],
        )

    def test_analysis_window_contains_full_cluster_and_flanks(self):
        self.assertEqual(locus_window(500_000, 900_000, 100_000),
                         (400_000, 1_000_000))
        self.assertEqual(locus_window(50_000, 300_000, 100_000),
                         (1, 400_000))
        with self.assertRaisesRegex(ValueError, "positive and ordered"):
            locus_window(300, 200, 100)

    def test_zero_locus_limit_selects_every_mapped_cluster(self):
        loci = [
            {"CHROM": "1", "start": 100, "end": 200, "lead_POS": 150,
             "lead_ID": "rs1", "liftover_status": "mapped"},
            {"CHROM": "1", "start": 500, "end": 600, "lead_POS": 550,
             "lead_ID": "rs2", "liftover_status": "excluded"},
            {"CHROM": "2", "start": 700, "end": 800, "lead_POS": 750,
             "lead_ID": "rs3", "liftover_status": "mapped"},
        ]
        significant = [
            {"CHROM": "1", "POS": 150, "ID": "rs1", "LP": 12},
            {"CHROM": "2", "POS": 750, "ID": "rs3", "LP": 10},
        ]

        selected = select_mapped_loci(loci, significant, max_loci=0)

        self.assertEqual([row["ID"] for row in selected], ["rs1", "rs3"])
        self.assertEqual(
            [(row["cluster_start"], row["cluster_end"]) for row in selected],
            [(100, 200), (700, 800)],
        )
        self.assertEqual(
            len(select_mapped_loci(loci, significant, max_loci=1)), 1
        )

    def test_liftover_maps_reverse_complemented_indel_span(self):
        lifter = StubLifter({
            ("chr2", 99): [("chr2", 200, "-", 7)],
            ("chr2", 100): [("chr2", 199, "-", 7)],
        })
        record = {"CHROM": "2", "POS": 100, "REF": "AC", "ALT": "A"}

        mapped, reason = map_variant(lifter, record)

        self.assertIsNone(reason)
        self.assertEqual(mapped["LD_POS"], 200)
        self.assertEqual(mapped["LD_REF"], "GT")
        self.assertEqual(mapped["LD_ALT"], "T")
        self.assertEqual(mapped["LIFTOVER_STRAND"], "-")

    def test_liftover_records_ambiguous_and_unmapped_variants(self):
        lifter = StubLifter({
            ("chr1", 99): [
                ("chr1", 199, "+", 1),
                ("chr1", 299, "+", 2),
            ]
        })
        ambiguous, reason = map_variant(
            lifter, {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G"}
        )
        unmapped, unmapped_reason = map_variant(
            lifter, {"CHROM": "1", "POS": 101, "REF": "A", "ALT": "G"}
        )

        self.assertIsNone(ambiguous)
        self.assertEqual(reason, "ambiguous_chain_mapping")
        self.assertIsNone(unmapped)
        self.assertEqual(unmapped_reason, "no_chain_mapping")

    def test_lifted_vcf_declares_target_contigs(self):
        lifter = StubLifter({("chr1", 99): [("chr1", 199, "+", 1)]})
        record = {
            "CHROM": "1", "POS": 100, "ID": "rs1", "REF": "A", "ALT": "G",
            "ES": 0.1, "SE": 0.05, "LP": 4.0, "SS": 1000,
        }
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "lifted.vcf"

            lifted, excluded = write_lifted_vcf([record], lifter, output)

            self.assertEqual(len(lifted), 1)
            self.assertEqual(excluded, [])
            self.assertIn("##contig=<ID=1>\n", output.read_text(encoding="utf-8"))

    def test_panel_match_collisions_keep_strongest_and_record_exclusions(self):
        normalized = [
            {"source_variant_key": "v1", "LD_CHROM": "1", "LD_POS": 200,
             "LD_REF": "A", "LD_ALT": "G", "LP": 4.8},
            {"source_variant_key": "v2", "LD_CHROM": "1", "LD_POS": 200,
             "LD_REF": "G", "LD_ALT": "A", "LP": 0.3},
            {"source_variant_key": "v3", "LD_CHROM": "1", "LD_POS": 201,
             "LD_REF": "C", "LD_ALT": "T", "LP": 2.0},
        ]
        matches = pd.DataFrame({
            "idx": [123, 123, -1],
            "flip": [False, True, False],
        })

        retained, excluded = resolve_panel_matches(normalized, matches)

        self.assertEqual([row["source_variant_key"] for row in retained], ["v1"])
        self.assertEqual(
            {(row["source_variant_key"], row["reason"]) for row in excluded},
            {("v2", "duplicate_panel_variant"), ("v3", "no_exact_panel_allele_match")},
        )
        self.assertEqual(
            next(row for row in excluded if row["source_variant_key"] == "v2")[
                "duplicate_of"
            ],
            "v1",
        )

    def test_matrix_export_normalizes_small_diagonal_drift_in_tsv_and_npy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            npz = root / "ld.npz"
            variants = root / "ld.variants.tsv"
            output_tsv = root / "matrix.tsv.gz"
            output_npy = root / "matrix.npy"
            np.savez(
                npz,
                matrix=np.array([[0.9999, 0.2], [0.2, 1.0001]]),
                idx=np.array([100, 101]),
            )
            variants.write_text(
                "idx\tchrom\tpos\tref\talt\trsid\tflipped\n"
                "100\t1\t200\tA\tG\trs1\tFalse\n"
                "101\t1\t201\tC\tT\trs2\tFalse\n",
                encoding="utf-8",
            )

            report = matrix_to_tsv(npz, variants, output_tsv)
            report_npy = matrix_to_npy(npz, variants, output_npy)

            self.assertEqual(report["variants"], 2)
            self.assertEqual(report_npy["variants"], 2)
            arr = np.load(output_npy)
            self.assertEqual(arr.shape, (2, 2))
            self.assertAlmostEqual(float(arr[0, 0]), 1.0)
            self.assertEqual(
                output_npy.with_suffix(".ids.txt").read_text(encoding="utf-8").splitlines(),
                ["1:200:A:G", "1:201:C:T"],
            )
            with gzip.open(output_tsv, "rt", encoding="utf-8") as stream:
                rows = list(csv.reader(stream, delimiter="\t"))
            self.assertEqual(rows[0][0], "LD_VARIANT_ID")
            self.assertEqual(rows[1][0], "1:200:A:G")

    def test_matrix_export_rejects_unavailable_ld(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            npz = root / "ld.npz"
            variants = root / "ld.variants.tsv"
            np.savez(npz, matrix=np.array([[1.0, np.nan], [np.nan, 1.0]]),
                     idx=np.array([100, 101]))
            variants.write_text(
                "idx\tchrom\tpos\tref\talt\trsid\tflipped\n"
                "100\t1\t200\tA\tG\trs1\tFalse\n"
                "101\t1\t201\tC\tT\trs2\tFalse\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "missing values"):
                matrix_to_npy(npz, variants, root / "matrix.npy")

    def test_region_method_runner_consumes_binary_npy_ld_and_names_coloc_matrices(self):
        project_root = Path(__file__).resolve().parents[2]
        environment = r_runtime_environment(project_root)
        rscript = shutil.which("Rscript", path=environment.get("PATH"))
        if not rscript:
            self.skipTest("Rscript is required")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            variants = root / "variants.tsv"
            matrix_npy = root / "matrix.npy"
            output = root / "gwas"
            cis_methods = output / "cis" / "methods"
            mhc_methods = output / "mhc" / "methods"
            records = [
                {
                    "CHROM": "1", "POS": 100 + i, "REF": ref, "ALT": alt,
                    "GWAS_VARIANT_ID": f"1:{1100 + i}:{ref}:{alt}",
                    "EZ": z, "LP": 8.0, "SS": 700 + i,
                }
                for i, (ref, alt, z) in enumerate([
                    ("A", "G", 5.5), ("C", "T", 5.2), ("G", "A", 0.3),
                    ("T", "C", -0.2), ("A", "C", 0.1),
                ])
            ]
            with variants.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=records[0].keys(),
                                        delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(records)

            ids = [f"1:{row['POS']}:{row['REF']}:{row['ALT']}" for row in records]
            ld = np.full((len(ids), len(ids)), 0.1, dtype="<f8")
            np.fill_diagonal(ld, 1.0)
            np.save(matrix_npy, ld, allow_pickle=False)
            matrix_npy.with_suffix(".ids.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")

            runner = project_root / "scripts" / "run_region_methods.R"
            for out_m, r_id, r_type in (
                (cis_methods, "cis_GENEA_chr1_100_200", "cis_gene"),
                (mhc_methods, "chr6_25500000_34000000_mhc", "mhc"),
            ):
                result = subprocess.run(
                    [
                        rscript, str(runner), "--variants", str(variants),
                        "--matrix", str(matrix_npy), "--outdir", str(out_m),
                        "--methods", "p_ld,susie_rss,susie_rrss",
                        "--reference-sample-size", "500",
                        "--region-id", r_id,
                        "--region-type", r_type,
                        "--write-tsv", "false",
                    ],
                    cwd=project_root,
                    env=environment,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)

            aggregator = runner.with_name("aggregate_method_results.R")
            result = subprocess.run(
                [rscript, str(aggregator), "--output-root", str(output)],
                cwd=project_root,
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            for filename in ("p-ld.rda", "susie-rss.rda", "susie-rrss.rda"):
                self.assertTrue((output / "methods" / filename).is_file())
            self.assertFalse((output / "methods" / "p-value.rda").exists())
            self.assertFalse((cis_methods / "susie_rss" / "variants.tsv").exists())

            check_r = (
                "load(commandArgs(TRUE)[1]); "
                "fit <- susie_rss$regions[['cis_GENEA_chr1_100_200']]$result$fit; "
                "stopifnot("
                "identical(susie_rss$method_id, 'susie_rss'), "
                "setequal(names(susie_rss$regions), c('cis_GENEA_chr1_100_200', 'chr6_25500000_34000000_mhc')), "
                "identical(colnames(fit$alpha)[1], '1:1100:A:G'), "
                "identical(colnames(fit$lbf_variable)[1], '1:1100:A:G'), "
                "identical(names(fit$pip)[1], '1:1100:A:G'))"
            )
            result = subprocess.run(
                [rscript, "-e", check_r, str(output / "methods" / "susie-rss.rda")],
                cwd=project_root,
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)


class InstalledPopulationIndexTests(unittest.TestCase):
    def setUp(self):
        import pyarrow as pa
        import pyarrow.parquet as pq

        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.panukb = self.tmp / "panukb"
        self.panukb.mkdir()
        table = pa.table({
            "contig": ["1", "1", "1", "2", "2"], "position": [100, 200, 900, 150, 250],
            "ref": ["A", "C", "G", "T", "A"], "alt": ["G", "T", "A", "C", "G"], "idx": [0, 1, 2, 3, 4],
        })
        self.whole = self.panukb / "UKBB.EUR.variant_index.parquet"
        pq.write_table(table, self.whole)
        Path(f"{self.whole}.json").write_text(json.dumps(
            {"population": "EUR", "region": "genome", "build": "GRCh37", "rows": 5}))

    def test_cut_keeps_only_the_region_and_the_columns(self):
        import pyarrow.parquet as pq

        out = self.tmp / "cut.parquet"
        self.assertEqual(cut_regional_index(self.whole, "1", 150, 300, out), 1)
        table = pq.read_table(out)
        self.assertEqual(table.column_names, ["contig", "position", "ref", "alt", "idx"])
        self.assertEqual(table.column("position").to_pylist(), [200])

    def test_regional_index_is_cut_from_the_installed_index_and_reused(self):
        import pyarrow.parquet as pq

        cache = self.tmp / "cache"
        index = ensure_regional_index("EUR", "1", 50, 300, cache, self.tmp / "no_scripts", self.panukb)
        self.assertEqual(pq.read_table(index).column("position").to_pylist(), [100, 200])
        self.assertEqual(json.loads(Path(f"{index}.json").read_text()),
                         {"population": "EUR", "region": "1:50-300", "build": "GRCh37"})
        self.assertEqual(ensure_regional_index("EUR", "1", 50, 300, cache, self.tmp / "no_scripts", self.panukb), index)

    def test_region_without_variants_is_an_error_and_leaves_no_index(self):
        cache = self.tmp / "cache"
        with self.assertRaisesRegex(ValueError, "no variants"):
            ensure_regional_index("EUR", "1", 400, 500, cache, self.tmp / "no_scripts", self.panukb)
        self.assertEqual(list((cache / "panukb").glob("*.parquet")), [])

    def test_missing_or_other_index_is_not_used(self):
        self.assertIsNone(installed_population_index("AFR", self.panukb))
        self.assertIsNone(installed_population_index("EUR", None))
        Path(f"{self.whole}.json").write_text(json.dumps(
            {"population": "EUR", "region": "1:1-10", "build": "GRCh37"}))
        with self.assertRaisesRegex(ValueError, "does not describe"):
            installed_population_index("EUR", self.panukb)


if __name__ == "__main__":
    unittest.main()
