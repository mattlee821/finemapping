"""The reference-data installer (scripts/reference_data.py): FASTA indexing, chain and FASTA checks, the gene annotation
(read back by run_finemapping's own reader), downloads, the Pan-UKB index record, and `check`. No network is used."""

import contextlib
import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import reference_data as rd
from run_finemapping import read_gene_annotation


def write_fasta(path: Path, records: dict[str, str], width: int = 60) -> Path:
    lines = []
    for name, sequence in records.items():
        lines.append(f">{name} description")
        lines += [sequence[i:i + width] for i in range(0, len(sequence), width)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def mapping_rows():
    columns = ["UNIPROT", "uniprot_gn_symbol", "ensembl_gene_id", "CHR", "START_hg38", "END_hg38"]
    rows = [
        ("P00001", "GENEA", "ENSG1", "1", "1000", "2000"),
        ("P00001", "GENEA", "", "1", "900", "1500"),          # same gene again, no Ensembl ID: merged into 1 span
        ("P00002", "GENEB", "ENSG2", "X", "5000", "6000"),
        ("P00003", "GENEC", "ENSG3", "2", "100", "200"),      # 1 accession, 2 genes: ambiguous
        ("P00003", "GENED", "ENSG4", "2", "300", "400"),
        ("P00004", "GENEE", "ENSG5", "", "", ""),             # no GRCh38 position: left out
        ("P00005", "GENEF", "ENSG6", "HSCHR6_MHC", "10", "20"),  # alternative contig: left out
    ]
    return [dict(zip(columns, row)) for row in rows]


class FastaTest(unittest.TestCase):
    def test_fai_matches_samtools_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = write_fasta(Path(tmp) / "a.fa", {"1": "ACGT" * 40, "2": "ACGT" * 10, "3": "A" * 60})
            self.assertEqual(rd.build_fai(fasta), 3)
            rows = [line.split("\t") for line in Path(f"{fasta}.fai").read_text().splitlines()]
            self.assertEqual(rows[0], ["1", "160", str(len(">1 description\n")), "60", "61"])
            self.assertEqual([r[:2] for r in rows], [["1", "160"], ["2", "40"], ["3", "60"]])
            offset_2 = int(rows[1][2])
            self.assertEqual(fasta.read_bytes()[offset_2:offset_2 + 4], b"ACGT")  # the offset is the first base of record 2

    def test_uneven_lines_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fa"
            path.write_text(">1\nACGT\nAC\nACGT\n", encoding="utf-8")
            with self.assertRaises(rd.ReferenceError):
                rd.build_fai(path)

    def test_check_fasta_wants_grch37(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "GRCh37.fa"
            fasta.write_text(">1\nACGT\n", encoding="utf-8")
            Path(f"{fasta}.fai").write_text("\n".join(f"{c}\t{rd.GRCH37_CHR1_LENGTH if c == '1' else 100}\t0\t60\t61"
                                                      for c in sorted(rd.STANDARD_CONTIGS)) + "\n", encoding="utf-8")
            with self.assertRaises(rd.ReferenceError) as raised:  # all contigs and the right chr1, but far too short in all
                rd.check_fasta(fasta)
            self.assertIn("truncated", str(raised.exception))
            Path(f"{fasta}.fai").write_text("1\t248956422\t0\t60\t61\n", encoding="utf-8")  # a GRCh38 chromosome 1
            with self.assertRaises(rd.ReferenceError):
                rd.check_fasta(fasta)


class ChainTest(unittest.TestCase):
    def test_chain_must_run_from_grch38_to_grch37(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "good.chain.gz"
            with gzip.open(good, "wt") as out:
                out.write(f"chain 1000 chr1 {rd.GRCH38_CHR1_LENGTH} + 0 100 chr1 {rd.GRCH37_CHR1_LENGTH} + 0 100 1\n100\n\n")
            rd.check_chain(good)
            reverse = Path(tmp) / "reverse.chain.gz"
            with gzip.open(reverse, "wt") as out:  # hg19ToHg38: the wrong way round
                out.write(f"chain 1000 chr1 {rd.GRCH37_CHR1_LENGTH} + 0 100 chr1 {rd.GRCH38_CHR1_LENGTH} + 0 100 1\n100\n\n")
            with self.assertRaises(rd.ReferenceError):
                rd.check_chain(reverse)
            garbage = Path(tmp) / "garbage.chain.gz"
            garbage.write_bytes(b"not a gzip file")
            with self.assertRaises(rd.ReferenceError):
                rd.check_chain(garbage)


class AnnotationTest(unittest.TestCase):
    def test_one_region_per_accession_and_readable_by_the_runner(self):
        annotation, ambiguous = rd.build_gene_annotation(mapping_rows(), "mapping test")
        self.assertEqual([a["uniprot_id"] for a in annotation], ["P00001", "P00002"])
        first = annotation[0]
        self.assertEqual((first["start"], first["end"], first["ensembl_id"], first["source_rows"]), (900, 2000, "ENSG1", 2))
        self.assertEqual(ambiguous, [{"uniprot_id": "P00003", "genes": "GENEC@2;GENED@2"}])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gene_regions.tsv"
            rd.write_tsv(path, annotation, rd.ANNOTATION_COLUMNS)
            gene = read_gene_annotation(path, "P00002")   # the runner's own reader
            self.assertEqual((gene["chrom"], gene["start"], gene["end"], gene["gene_symbol"]), ("X", 5000, 6000, "GENEB"))
            with self.assertRaises(ValueError):           # an ambiguous accession is absent, so the runner stops there
                read_gene_annotation(path, "P00003")

    def test_install_refuses_a_mapping_that_is_too_small(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "m.tsv"
            pd.DataFrame(mapping_rows()).to_csv(table, sep="\t", index=False)
            with mock.patch.object(rd, "export_usefuldata", lambda rda, destination, rscript: destination.write_text(table.read_text())):
                with self.assertRaises(rd.ReferenceError):
                    rd.install_annotation(Path(tmp) / "root", force=True, rda=table)


class DownloadTest(unittest.TestCase):
    def test_fetch_a_local_file_and_fail_on_a_missing_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.txt"
            source.write_text("content", encoding="utf-8")
            destination = rd.fetch(source.as_uri(), Path(tmp) / "sub" / "copy.txt", retries=1)
            self.assertEqual(destination.read_text(), "content")
            self.assertFalse(list((Path(tmp) / "sub").glob("*.part")))
            with mock.patch.object(rd.time, "sleep"), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(rd.ReferenceError):
                    rd.fetch((Path(tmp) / "absent.txt").as_uri(), Path(tmp) / "never.txt", retries=2)
            self.assertFalse((Path(tmp) / "never.txt").exists())


class PanukbTest(unittest.TestCase):
    def write_index(self, root: Path, population: str, rows: int, recorded: int | None = None) -> Path:
        index, sidecar = rd.panukb_paths(root, population)
        index.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"contig": ["1"] * rows, "position": range(1, rows + 1), "ref": ["A"] * rows, "alt": ["G"] * rows,
                      "idx": range(rows)}).to_parquet(index)
        sidecar.write_text(json.dumps({"population": population, "region": "genome", "build": "GRCh37",
                                       "rows": rows if recorded is None else recorded}), encoding="utf-8")
        return index

    def test_check_reads_the_index_against_its_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_index(root, "EUR", 5)
            self.assertEqual(rd.check_panukb(root, "EUR"), 5)
            self.write_index(root, "AFR", 5, recorded=9)
            with self.assertRaises(rd.ReferenceError):
                rd.check_panukb(root, "AFR")
            with self.assertRaises(rd.ReferenceError):
                rd.check_panukb(root, "EAS")  # not installed

    def test_install_builds_missing_populations_only_and_clears_stale_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_index(root, "EUR", 5)  # already installed: not rebuilt
            stale = Path(f"{rd.panukb_paths(root, 'AFR')[0]}.export.tsv.bgz")
            stale.parent.mkdir(parents=True, exist_ok=True)
            stale.write_text("interrupted")
            calls = []

            def fake_build(command, check, env):
                calls.append(command[2:])
                self.assertFalse(stale.exists())  # the builder refuses to run beside an earlier export
                self.assertEqual(env["GWASFM_HAIL_DRIVER_MEMORY"], "8g")
                index = Path(command[3])
                self.write_index(root, command[2], 7)
                index.touch()
                return mock.Mock(returncode=0)

            with mock.patch.object(rd.subprocess, "run", fake_build), contextlib.redirect_stdout(io.StringIO()):
                rd.install_panukb(root, ("EUR", "AFR"), force=False)
            self.assertEqual([c[0] for c in calls], ["AFR"])
            self.assertEqual(rd.check_panukb(root, "AFR"), 7)
            self.assertEqual(json.loads(rd.panukb_paths(root, "AFR")[1].read_text())["region"], "genome")
            self.assertIn("panukb/UKBB.AFR.variant_index.parquet", rd.read_manifest(root)["items"])


class CheckCommandTest(unittest.TestCase):
    def test_check_reports_what_is_missing_and_what_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = mock.Mock(root=str(root), only="annotation", populations=None)
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(rd.command_check(args), 1)
            self.assertIn("MISSING", out.getvalue())
            annotation, _ = rd.build_gene_annotation(mapping_rows(), "mapping test")
            path = root / "annotation" / "gene_regions_GRCh38.tsv"
            rd.write_tsv(path, annotation, rd.ANNOTATION_COLUMNS)
            rd.record(root, "annotation/gene_regions_GRCh38.tsv", source="test")
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(rd.command_check(args), 0)
            path.write_text(path.read_text() + "\n")
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(rd.command_check(args), 1)
            self.assertIn("CHANGED", out.getvalue())

    def test_paths_for_nextflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = mock.Mock(root=tmp, nextflow=True)
            with contextlib.redirect_stdout(io.StringIO()) as out:
                rd.command_paths(args)
            self.assertIn(f"--chain {Path(tmp).resolve()}/genome/hg38ToHg19.over.chain.gz", out.getvalue())
            self.assertIn("--gene_annotation", out.getvalue())

    def test_environment_root(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {"GWASFM_REFERENCE_ROOT": tmp}):
            self.assertEqual(rd.reference_root(None), Path(tmp).resolve())


class RscriptTest(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_export_reads_an_rda_with_real_r(self):
        with tempfile.TemporaryDirectory() as tmp:
            rda = Path(tmp) / "mapping.rda"
            code = (f"{rd.USEFULDATA_DATASET} <- data.frame(UNIPROT = c('P1', 'P2'), uniprot_gn_symbol = c('G1', 'G2'), "
                    f"ensembl_gene_id = c('E1', NA), CHR = c('1', 'X'), START_hg38 = c(10, 20), END_hg38 = c(30, 40)); "
                    f"save({rd.USEFULDATA_DATASET}, file = '{rda}')")
            subprocess.run(["Rscript", "--vanilla", "-e", code], check=True, capture_output=True)
            exported = Path(tmp) / "mapping.tsv"
            rd.export_usefuldata(rda, exported, "Rscript")
            with exported.open(encoding="utf-8") as stream:
                table = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([r["UNIPROT"] for r in table], ["P1", "P2"])
            annotation, _ = rd.build_gene_annotation(table, "test")
            self.assertEqual([(a["chrom"], a["start"], a["end"], a["ensembl_id"]) for a in annotation], [("1", 10, 30, "E1"), ("X", 20, 40, "")])
            with self.assertRaises(rd.ReferenceError):
                rd.export_usefuldata(Path(tmp) / "absent.rda", Path(tmp) / "x.tsv", "Rscript")


if __name__ == "__main__":
    unittest.main()
