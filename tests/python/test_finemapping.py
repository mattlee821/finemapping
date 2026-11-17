import csv
import gzip
import os
import shutil
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from run_finemapping import (
    cluster_significant,
    cut_regional_index,
    ensure_regional_index,
    installed_population_index,
    locus_window,
    map_variant,
    matrix_to_tsv,
    resolve_panel_matches,
    select_mapped_loci,
    write_lifted_vcf,
)


class StubLifter:
    def __init__(self, mapping):
        self.mapping = mapping

    def convert_coordinate(self, chromosome, position):
        return self.mapping.get((chromosome, position), [])


class FineMappingTests(unittest.TestCase):
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

    def test_matrix_export_normalizes_small_diagonal_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            npz = root / "ld.npz"
            variants = root / "ld.variants.tsv"
            output = root / "matrix.tsv.gz"
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

            report = matrix_to_tsv(npz, variants, output)

            self.assertEqual(report["variants"], 2)
            with gzip.open(output, "rt", encoding="utf-8") as stream:
                rows = list(csv.reader(stream, delimiter="\t"))
            self.assertEqual(rows[0][0], "LD_VARIANT_ID")
            self.assertEqual(rows[1][0], "1:200:A:G")
            self.assertAlmostEqual(float(rows[1][1]), 1.0)
            self.assertAlmostEqual(float(rows[2][2]), 1.0)

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
                matrix_to_tsv(npz, variants, root / "matrix.tsv.gz")

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_region_method_runner_consumes_aligned_ld_matrix(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            variants = root / "variants.tsv"
            matrix = root / "matrix.tsv.gz"
            output = root / "gwas"
            gene_methods = output / "gene_region" / "methods"
            locus_methods = output / "loci" / "chr1_100_200" / "methods"
            records = [
                {"CHROM": "1", "POS": 100 + i, "REF": ref, "ALT": alt,
                 "EZ": z, "LP": 2, "SS": 700 + i}
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
            ld = np.full((len(ids), len(ids)), 0.1)
            np.fill_diagonal(ld, 1)
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(["LD_VARIANT_ID", *ids])
                for variant_id, row in zip(ids, ld):
                    writer.writerow([variant_id, *row])

            runner = Path(__file__).resolve().parents[2] / "scripts" / "run_region_methods.R"
            environment = os.environ.copy()
            tools_root = Path(
                environment.get(
                    "GWASFM_TOOLS_ROOT",
                    str(Path.home() / ".cache" / "gwasfinemapping"),
                )
            ).expanduser()
            environment["RENV_PATHS_LIBRARY"] = str(tools_root / "renv-library")
            result = subprocess.run(
                [
                    "Rscript", str(runner), "--variants", str(variants),
                    "--matrix", str(matrix), "--outdir", str(gene_methods),
                    "--methods", "p_value,p_ld,susie_rss,susie_rrss",
                    "--reference-sample-size", "500",
                    "--region-id", "gene_region",
                    "--region-type", "gene_region",
                ],
                cwd=runner.parents[1],
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            result = subprocess.run(
                [
                    "Rscript", str(runner), "--variants", str(variants),
                    "--matrix", str(matrix), "--outdir", str(locus_methods),
                    "--methods", "p_value,p_ld,susie_rss,susie_rrss",
                    "--reference-sample-size", "500",
                    "--region-id", "chr1_100_200",
                    "--region-type", "genome_wide_locus",
                ],
                cwd=runner.parents[1],
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            aggregator = runner.with_name("aggregate_method_results.R")
            result = subprocess.run(
                [
                    "Rscript", str(aggregator), "--output-root", str(output),
                ],
                cwd=runner.parents[1],
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            self.assertTrue((gene_methods / "method_summary.tsv").is_file())
            self.assertTrue((gene_methods / "susie_rss" / "variants.tsv").is_file())
            self.assertTrue(
                (gene_methods / "susie_rrss" / "credible_sets.tsv").is_file()
            )
            for filename in (
                "p-value.rda", "p-ld.rda", "susie-rss.rda", "susie-rrss.rda"
            ):
                self.assertTrue((output / "methods" / filename).is_file())
            self.assertFalse(
                (gene_methods / "susie_rss" / "region_result.rds").exists()
            )
            result = subprocess.run(
                [
                    "Rscript", "-e",
                    "load(commandArgs(TRUE)[1]); "
                    "stopifnot(identical(susie_rss$method_id, 'susie_rss'), "
                    "setequal(names(susie_rss$regions), "
                    "c('gene_region', 'chr1_100_200')), "
                    "identical(susie_rss$regions$gene_region$region_type, "
                    "'gene_region'), "
                    "identical(susie_rss$regions$chr1_100_200$region_type, "
                    "'genome_wide_locus'), "
                    "is.list(susie_rss$regions$gene_region$result$fit))",
                    str(output / "methods" / "susie-rss.rda"),
                ],
                cwd=runner.parents[1],
                env=environment,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)


class InstalledPopulationIndexTests(unittest.TestCase):
    def setUp(self):
        import json
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
        import json
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
        import json

        self.assertIsNone(installed_population_index("AFR", self.panukb))
        self.assertIsNone(installed_population_index("EUR", None))
        Path(f"{self.whole}.json").write_text(json.dumps(
            {"population": "EUR", "region": "1:1-10", "build": "GRCh37"}))
        with self.assertRaisesRegex(ValueError, "does not describe"):
            installed_population_index("EUR", self.panukb)


if __name__ == "__main__":
    unittest.main()
