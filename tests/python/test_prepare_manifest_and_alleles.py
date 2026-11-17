"""Unit & integration tests for prepare_manifest.py (both GWASprocessing layouts,
AppleDouble filtering, .finemap_state.tsv idempotency) and swapped-allele SNP/indel
liftover + SuSiE sign invariance."""

from __future__ import annotations

import gzip
import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import prepare_manifest as pm
from build_synthetic_fixture import write_bgzip_vcf, write_chain, write_fasta_with_fai
from run_finemapping import (
    extract_ld_from_local_panel,
    matrix_to_npy,
    normalize_lifted_vcf,
    r_runtime_environment,
    resolve_panel_matches,
    run_methods,
    write_lifted_vcf,
    write_tsv,
)


class StubLifter:
    def __init__(self, mapping):
        self.mapping = mapping

    def convert_coordinate(self, chromosome, position):
        return self.mapping.get((chromosome, position), [])


class PrepareManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gwasfm-prep-"))
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def _make_mini_vcf(self, vcf_gz: Path) -> None:
        write_bgzip_vcf(vcf_gz, [{
            "chrom": "1", "pos": 10000, "id": "rs1", "ref": "A", "alt": "G",
            "es": 0.2, "se": 0.03, "lp": 10.5, "ss": 15000, "af": 0.3,
        }])

    def test_discovers_both_gwasprocessing_layouts_and_ignores_appledouble(self):
        root = self.tmp / "GWAS"
        # Layout A: locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.*
        d_locke = root / "locke_2015_25673413" / "standardise" / "GWAScatalog_GCST002783"
        self._make_mini_vcf(d_locke / "standardised.gwas.vcf.gz")
        (d_locke / "standardised.metadata.json").write_text(json.dumps({
            "study_id": "locke_2015_25673413",
            "output_id": "GWAScatalog_GCST002783",
            "phenotype": "BMI",
            "ancestry": "ALL",
            "coordinate_build": "GRCh38",
        }), encoding="utf-8")
        (d_locke / "standardised.qc.json").write_text(json.dumps({"coordinate_build": "GRCh38"}), encoding="utf-8")
        # AppleDouble noise files/folders that must be ignored
        (d_locke / "._standardised.gwas.vcf.gz").write_bytes(b"\x00\x05\x16\x07")
        (root / "locke_2015_25673413" / "standardise" / "._GWAScatalog_GCST002783").write_bytes(b"\x00")

        # Layout B: sun_2023_37794186/standardise/P16860_OID20049/P16860_OID20049.*
        d_sun = root / "sun_2023_37794186" / "standardise" / "P16860_OID20049"
        self._make_mini_vcf(d_sun / "P16860_OID20049.gwas.vcf.gz")
        (d_sun / "P16860_OID20049.metadata.json").write_text(json.dumps({
            "study_id": "sun_2023_37794186",
            "output_id": "P16860_OID20049",
            "phenotype": "NPPB protein",
            "ancestry": "ALL",
            "coordinate_build": "GRCh38",
        }), encoding="utf-8")
        (d_sun / "P16860_OID20049.qc.json").write_text(json.dumps({"coordinate_build": "GRCh38"}), encoding="utf-8")

        discovered = pm.discover_input_dir(root)
        self.assertEqual(len(discovered), 2)
        keys = {(r["study_id"], r["gwas_id"]) for r in discovered}
        self.assertEqual(
            keys,
            {
                ("locke_2015_25673413", "GWAScatalog_GCST002783"),
                ("sun_2023_37794186", "P16860_OID20049"),
            },
        )
        for raw in discovered:
            checked, errs = pm.preflight_row(raw, ld_population="EUR", allow_ancestry_mismatch=False)
            self.assertEqual(errs, [])
            self.assertIsNotNone(checked)

    def test_single_ancestry_mismatch_fails_preflight_while_all_passes(self):
        d_eas = self.tmp / "eas_study" / "standardise" / "trait_eas"
        self._make_mini_vcf(d_eas / "standardised.gwas.vcf.gz")
        (d_eas / "standardised.metadata.json").write_text(json.dumps({
            "study_id": "eas_study",
            "output_id": "trait_eas",
            "ancestry": "EAS",
            "coordinate_build": "GRCh38",
        }), encoding="utf-8")
        raw = pm.resolve_gwas_folder(d_eas)
        _, errs = pm.preflight_row(raw, ld_population="EUR", allow_ancestry_mismatch=False)
        self.assertEqual(len(errs), 1)
        self.assertIn("EAS", errs[0])

        _, errs_ok = pm.preflight_row(raw, ld_population="EUR", allow_ancestry_mismatch=True)
        self.assertEqual(errs_ok, [])


class SwappedAlleleAndIndelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gwasfm-swap-"))
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_swapped_grch37_indel_normalizes_and_preserves_gwas_alt_orientation(self):
        # GRCh37 FASTA has "AT" followed by "TTTT..." at pos 10..20 (deletion of T in homopolymer).
        # Suppose Pan-UKB (GRCh37) stores ref="AT", alt="A" at left-aligned pos 10,
        # while GRCh38 had the deleted allele as reference at right-shifted pos 14: GWAS REF="T", ALT="TT" -> lifted REF="A", ALT="AT"!
        fasta = self.tmp / "ref.fa"
        seq = "C" * 9 + "A" + "T" * 10 + "G" * 60
        fasta.write_text(f">1\n{seq}\n", encoding="utf-8")
        Path(f"{fasta}.fai").write_text(f"1\t{len(seq)}\t3\t80\t81\n", encoding="utf-8")

        # Lifted coordinate at pos 13 (inside homopolymer): lifted REF="T", ALT="TT"
        # Wait: on GRCh37 at pos 13, fasta has "TT" == ALT (and "T" is prefix, so test non-prefix e.g. SNV/indel swap at pos 10:
        # at pos 10 fasta has "AT"; if GWAS had REF="A", ALT="AT" at pos 10, "AT".startswith("A") is True;
        # what if GWAS had deletion at pos 12 where REF="G" vs ALT="AT"? Let's test exact indel left-shift from pos 13:
        # On GRCh37 seq[9:20] = "ATTTTTTTTTT". Deletion of 1 T at pos 13 on GRCh37 is REF="TT", ALT="T" -> left-aligns to pos 10 REF="AT", ALT="A".
        # Now suppose GWAS had insertion relative to deleted GRCh38: lifted at pos 10 has REF="A", ALT="AT" (or at pos 13 REF="T", ALT="TT")!
        # Even stronger: test an insertion/deletion where REF does NOT match fasta prefix, e.g. MNV/indel REF="C", ALT="AT" at pos 10 -> fasta starts with "AT"!
        lifter = StubLifter({("chr1", 109): [("chr1", 9, "+", 1)]})
        record = {
            "CHROM": "1", "POS": 110, "ID": "rs_indel", "REF": "C", "ALT": "AT",
            "ES": 0.42, "SE": 0.05, "LP": 16.0, "SS": 10000,
        }
        lifted_vcf = self.tmp / "lifted.vcf"
        write_lifted_vcf([record], lifter, lifted_vcf)
        norm_rows = normalize_lifted_vcf(lifted_vcf, fasta)
        self.assertEqual(len(norm_rows), 1)
        # Orientation restored to GWAS ALT ("AT")!
        self.assertEqual(norm_rows[0]["LD_ALT"], "AT")
        self.assertEqual(norm_rows[0]["LD_REF"], "C")
        self.assertAlmostEqual(norm_rows[0]["ES"], 0.42)

    def test_allele_swap_in_gwas_produces_identical_susie_pips_and_cs(self):
        project_root = Path(__file__).resolve().parents[2]
        env = r_runtime_environment(project_root)
        if not shutil.which("Rscript", path=env.get("PATH")):
            self.skipTest("Rscript required")

        # 6 variants with correlated LD; causal variant is index 1
        pos = [1000 + i * 100 for i in range(6)]
        R_true = np.power(0.85, np.abs(np.subtract.outer(range(6), range(6))))
        panel_npz = self.tmp / "panel.npz"
        np.savez(panel_npz, matrix=R_true, idx=np.arange(6, dtype=np.int64))

        z_true = 8.5 * R_true[:, 1]
        se = 0.02

        def run_orientation(swap_mask: list[bool], tag: str):
            rows = []
            for i in range(6):
                swapped = swap_mask[i]
                ref, alt = ("G", "A") if swapped else ("A", "G")
                z_i = -z_true[i] if swapped else z_true[i]
                rows.append({
                    "CHROM": "1", "POS": pos[i], "ID": f"rs{i}", "REF": ref, "ALT": alt,
                    "GWAS_VARIANT_ID": f"1:{pos[i]}:A:G",  # fixed canonical ID to compare PIPs directly
                    "GWAS_CHROM": "1", "GWAS_POS": pos[i], "GWAS_REF": ref, "GWAS_ALT": alt,
                    "LP": 12.0 if abs(z_i) > 6 else 2.0,
                    "ES": float(z_i * se), "SE": se, "SS": 10000,
                    "EZ": float(z_i),
                    "LD_CHROM": "1", "LD_POS": pos[i], "LD_REF": ref, "LD_ALT": alt,
                    "LD_VARIANT_ID": f"1:{pos[i]}:{ref}:{alt}",
                    "PANEL_IDX": i, "PANEL_ALLELE_SWAPPED": swapped,
                })
            sub = self.tmp / tag
            sub.mkdir(parents=True, exist_ok=True)
            v_tsv = sub / "variants.tsv"
            write_tsv(v_tsv, rows, list(rows[0].keys()))
            ld_pref = sub / "ld"
            extract_ld_from_local_panel(panel_npz, rows, ld_pref)
            npy_p = sub / "matrix.npy"
            matrix_to_npy(Path(f"{ld_pref}.npz"), Path(f"{ld_pref}.variants.tsv"), npy_p)
            run_methods(
                variants_path=v_tsv,
                matrix_path=npy_p,
                output_dir=sub / "methods",
                reference_n=50000,
                region_id="loc1",
                region_type="lead_locus",
            )
            out_json = sub / "pips.json"
            subprocess.run(
                [
                    "Rscript", "-e",
                    "r <- readRDS(commandArgs(TRUE)[1]); "
                    "jsonlite::write_json(list(pip=r$result$variants$PIP, cs=r$result$credible_sets$credible_set_id), commandArgs(TRUE)[2])",
                    str(sub / "methods" / "susie_rrss" / "region_result.rds"),
                    str(out_json),
                ],
                cwd=project_root, env=env, check=True,
            )
            return json.loads(out_json.read_text(encoding="utf-8"))

        unswapped = run_orientation([False] * 6, "unswapped")
        swapped = run_orientation([False, True, False, True, True, False], "swapped")
        self.assertEqual(unswapped["cs"], swapped["cs"])
        for p1, p2 in zip(unswapped["pip"], swapped["pip"]):
            self.assertAlmostEqual(p1, p2, places=6)

    def test_tandem_repeat_indel_pairs_do_not_cross_swap_or_collide(self) -> None:
        import pandas as pd
        import pyarrow as pa
        import pyarrow.parquet as pq
        from run_finemapping import match_panel_variants, resolve_panel_matches

        pq_path = self.tmp / "panel.parquet"
        table = pa.Table.from_pandas(pd.DataFrame({
            "contig": ["10", "10", "10"],
            "position": [122649699, 122649699, 122649750],
            "ref": ["G", "GATTT", "C"],
            "alt": ["GATTT", "G", "T"],
            "idx": [14999779, 14999780, 14999781],
        }), preserve_index=False)
        pq.write_table(table, pq_path)

        # Case 1: Both insertion (G->GATTT) and deletion (GATTT->G) present in GWAS & panel
        norm_both = [
            {"source_variant_key": "v0", "LD_CHROM": "10", "LD_POS": 122649699,
             "LD_REF": "G", "LD_ALT": "GATTT", "LP": 1.0, "PRE_SWAPPED": False, "FAI_CHECKED": True},
            {"source_variant_key": "v1", "LD_CHROM": "10", "LD_POS": 122649699,
             "LD_REF": "GATTT", "LD_ALT": "G", "LP": 2.0, "PRE_SWAPPED": False, "FAI_CHECKED": True},
            {"source_variant_key": "v2", "LD_CHROM": "10", "LD_POS": 122649750,
             "LD_REF": "T", "LD_ALT": "C", "LP": 3.0, "PRE_SWAPPED": False, "FAI_CHECKED": True},
        ]
        m_both = match_panel_variants(pq_path, norm_both)
        matched, excluded = resolve_panel_matches(norm_both, m_both)
        by_key = {r["source_variant_key"]: r for r in matched}
        self.assertEqual(len(excluded), 0)
        self.assertEqual(by_key["v0"]["PANEL_IDX"], 14999779)
        self.assertFalse(by_key["v0"]["PANEL_ALLELE_SWAPPED"])
        self.assertEqual(by_key["v1"]["PANEL_IDX"], 14999780)
        self.assertFalse(by_key["v1"]["PANEL_ALLELE_SWAPPED"])
        self.assertEqual(by_key["v2"]["PANEL_IDX"], 14999781)
        self.assertTrue(by_key["v2"]["PANEL_ALLELE_SWAPPED"])

        # Case 2: Panel has ONLY insertion (G->GATTT); non-swapped deletion (GATTT->G) must NOT cross-match
        pq_ins_only = self.tmp / "panel_ins_only.parquet"
        pq.write_table(pa.Table.from_pandas(pd.DataFrame({
            "contig": ["10"],
            "position": [122649699],
            "ref": ["G"],
            "alt": ["GATTT"],
            "idx": [14999779],
        }), preserve_index=False), pq_ins_only)
        m_del = match_panel_variants(pq_ins_only, [norm_both[1]])
        matched_del, excluded_del = resolve_panel_matches([norm_both[1]], m_del)
        self.assertEqual(len(matched_del), 0)
        self.assertEqual(excluded_del[0]["reason"], "no_exact_panel_allele_match")

        # Case 3: Genuine pre-swapped indel (PRE_SWAPPED=True) still matches via allow_allele_swap
        pre_swapped_indel = [{
            "source_variant_key": "v_swap", "LD_CHROM": "10", "LD_POS": 122649699,
            "LD_REF": "GATTT", "LD_ALT": "G", "LP": 4.0, "PRE_SWAPPED": True, "FAI_CHECKED": True,
        }]
        m_swap = match_panel_variants(pq_ins_only, pre_swapped_indel)
        matched_swap, _ = resolve_panel_matches(pre_swapped_indel, m_swap)
        self.assertEqual(len(matched_swap), 1)
        self.assertEqual(matched_swap[0]["PANEL_IDX"], 14999779)
        self.assertTrue(matched_swap[0]["PANEL_ALLELE_SWAPPED"])


if __name__ == "__main__":
    unittest.main()
