import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from build_panukb_variant_index import spark_conf, write_index_from_export


class VariantIndexBuilderTests(unittest.TestCase):
    def write_export(self, path: Path, rows: list[str]) -> None:
        with gzip.open(path, "wt", encoding="utf-8") as export:
            export.write("contig\tposition\tref\talt\tidx\tallele_count\n")
            export.writelines(f"{row}\n" for row in rows)

    def test_writes_export_in_chunks_with_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            export = Path(temp) / "variants.tsv.bgz"
            output = Path(temp) / "variants.parquet"
            self.write_export(
                export,
                [
                    "1\t100\tA\tC\t9\t2",
                    "1\t200\tG\tT\t3\t2",
                    "2\t300\tC\tA\t7\t2",
                ],
            )

            count = write_index_from_export(str(export), str(output), chunk_rows=2)
            table = pq.read_table(output)

            self.assertEqual(count, 3)
            self.assertEqual(table.column_names, ["contig", "position", "ref", "alt", "idx"])
            self.assertEqual(table.column("idx").to_pylist(), [9, 3, 7])

    def test_rejects_non_biallelic_rows_and_removes_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            export = Path(temp) / "variants.tsv.bgz"
            output = Path(temp) / "variants.parquet"
            self.write_export(
                export,
                [
                    "1\t100\tA\tC\t9\t2",
                    "1\t200\tG\tT\t3\t3",
                ],
            )

            with self.assertRaisesRegex(ValueError, "Non-biallelic"):
                write_index_from_export(str(export), str(output), chunk_rows=1)

            self.assertFalse(output.exists())
            self.assertFalse(Path(f"{output}.tmp").exists())


    def test_spark_conf_uses_the_loopback_address(self) -> None:
        # a node named cn17_mgmt gives Spark the invalid URL spark://HeartbeatReceiver@cn17_mgmt:46155 unless the driver
        # address is set
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("GWASFM_SPARK_DRIVER_HOST", None)
            os.environ.pop("GWASFM_HAIL_DRIVER_MEMORY", None)
            conf = spark_conf()
        self.assertEqual((conf["spark.driver.host"], conf["spark.driver.bindAddress"]), ("127.0.0.1", "127.0.0.1"))
        self.assertEqual(conf["spark.driver.memory"], "2g")
        self.assertIn("hadoop-aws", conf["spark.jars.packages"])
        with mock.patch.dict(os.environ, {"GWASFM_SPARK_DRIVER_HOST": "10.0.0.5", "GWASFM_HAIL_DRIVER_MEMORY": "8g"}):
            conf = spark_conf()
        self.assertEqual((conf["spark.driver.host"], conf["spark.driver.memory"]), ("10.0.0.5", "8g"))


if __name__ == "__main__":
    unittest.main()
