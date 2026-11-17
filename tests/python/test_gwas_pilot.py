import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from gwas_pilot import discover_inputs, preflight_study, write_manifest, write_summary


class GwasPilotTests(unittest.TestCase):
    def make_study(self, root: Path, study_id: str) -> dict[str, Path]:
        directory = root / study_id
        directory.mkdir()
        paths = {
            "vcf": directory / f"{study_id}.gwas.vcf.gz",
            "index": directory / f"{study_id}.gwas.vcf.gz.tbi",
            "metadata": directory / f"{study_id}.metadata.json",
            "qc": directory / f"{study_id}.qc.json",
        }
        with gzip.open(paths["vcf"], "wt", encoding="utf-8") as stream:
            stream.write(
                "##fileformat=VCFv4.2\n"
                "##reference=GRCh38.fa\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                f"{study_id}\n"
            )
        paths["index"].write_bytes(b"index")
        paths["metadata"].write_text(
            json.dumps(
                {
                    "output_id": study_id,
                    "standard_format": "GWAS-VCF",
                    "coordinate_build": "GRCh38",
                    "ancestry": "EUR",
                    "phenotype": "test trait",
                }
            ),
            encoding="utf-8",
        )
        paths["qc"].write_text(
            json.dumps(
                {
                    "coordinate_build": "GRCh38",
                    "build_detection": {"build": "hg38"},
                }
            ),
            encoding="utf-8",
        )
        return paths

    def test_manifest_is_sorted_capped_and_reports_discovery_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_study(root, "P43320")
            self.make_study(root, "P04049")
            output = root / "manifest.tsv"
            summary_path = root / "discovery.json"

            summary = write_manifest(str(root), str(output), str(summary_path), 1)

            with output.open(encoding="utf-8", newline="") as stream:
                entries = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([entry["study_id"] for entry in entries], ["P04049"])
            self.assertEqual(summary["discovered_count"], 2)
            self.assertEqual(summary["selected_count"], 1)
            self.assertEqual(json.loads(summary_path.read_text())["max_inputs"], 1)

    def test_manifest_fails_explicitly_for_missing_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.make_study(root, "P43320")
            paths["index"].unlink()

            with self.assertRaisesRegex(ValueError, "missing index"):
                discover_inputs(str(root), 1)

    def test_manifest_refuses_to_resume_with_different_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_study(root, "P43320")
            existing = root / "existing.tsv"
            existing.write_text("study_id\tvcf\nold\told.vcf\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "existing output manifest differs"):
                write_manifest(
                    str(root),
                    str(root / "manifest.tsv"),
                    str(root / "summary.json"),
                    1,
                    existing_manifest=str(existing),
                )

    def test_preflight_validates_build_and_records_index_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.make_study(root, "P43320")
            output = root / "report.json"
            with patch(
                "gwas_pilot.subprocess.run",
                return_value=CompletedProcess(
                    args=["bcftools", "index"],
                    returncode=0,
                    stdout="25359359\n",
                    stderr="",
                ),
            ) as run:
                report = preflight_study(
                    "P43320",
                    str(paths["vcf"]),
                    str(paths["index"]),
                    str(paths["metadata"]),
                    str(paths["qc"]),
                    str(output),
                    "/source/path/P43320.gwas.vcf.gz",
                    "/source/path/P43320.gwas.vcf.gz.tbi",
                    "/source/path/P43320.metadata.json",
                    "/source/path/P43320.qc.json",
                )

            run.assert_called_once()
            self.assertEqual(report["index_record_count"], 25_359_359)
            self.assertEqual(report["genome_build"], "GRCh38")
            self.assertEqual(report["vcf"], "/source/path/P43320.gwas.vcf.gz")
            self.assertEqual(json.loads(output.read_text())["sample_id"], "P43320")

    def test_preflight_rejects_reference_build_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.make_study(root, "P43320")
            qc = json.loads(paths["qc"].read_text())
            qc["build_detection"]["build"] = "hg19"
            paths["qc"].write_text(json.dumps(qc), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "consistently identify GRCh38"):
                preflight_study(
                    "P43320",
                    str(paths["vcf"]),
                    str(paths["index"]),
                    str(paths["metadata"]),
                    str(paths["qc"]),
                    str(root / "report.json"),
                )

    def test_preflight_rejects_invalid_build_detection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self.make_study(root, "P43320")
            qc = json.loads(paths["qc"].read_text())
            qc["build_detection"] = []
            paths["qc"].write_text(json.dumps(qc), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "QC build_detection is missing"):
                preflight_study(
                    "P43320",
                    str(paths["vcf"]),
                    str(paths["index"]),
                    str(paths["metadata"]),
                    str(paths["qc"]),
                    str(root / "report.json"),
                )

    def test_summary_is_sorted_and_has_a_header(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "P43320.json"
            second = root / "P04049.json"
            reports = []
            for study_id, report_path in (("P43320", first), ("P04049", second)):
                report_path.write_text(
                    json.dumps(
                        {
                            "study_id": study_id,
                            "ancestry": "EUR",
                            "genome_build": "GRCh38",
                            "sample_id": study_id,
                            "index_record_count": 10,
                            "vcf_bytes": 100,
                            "index_bytes": 20,
                            "phenotype": "test",
                        }
                    ),
                    encoding="utf-8",
                )
                reports.append(str(report_path))
            output = root / "summary.tsv"

            write_summary(reports, str(output))

            with output.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([row["study_id"] for row in rows], ["P04049", "P43320"])
            self.assertIn("index_record_count", rows[0])


if __name__ == "__main__":
    unittest.main()
