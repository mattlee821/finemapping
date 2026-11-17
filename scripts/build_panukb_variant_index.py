#!/usr/bin/env python3
"""Build a bounded-memory Parquet index for a public Pan-UKB LD panel."""

import argparse
import os
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ESSENTIAL_COLUMNS = ["contig", "position", "ref", "alt", "idx"]
EXPORT_COLUMNS = ESSENTIAL_COLUMNS + ["allele_count"]
CHUNK_ROWS = 250_000
HADOOP_AWS_VERSION = "3.3.4"
POPULATIONS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")


def write_index_from_export(export_path: str, output_path: str, chunk_rows: int = CHUNK_ROWS) -> int:
    temporary_output = f"{output_path}.tmp"
    if os.path.exists(output_path) or os.path.exists(temporary_output):
        raise FileExistsError(f"Refusing to overwrite an existing index: {output_path}")

    writer = None
    row_count = 0
    try:
        chunks = pd.read_csv(
            export_path,
            sep="\t",
            compression="gzip",
            dtype={
                "contig": str,
                "position": "int64",
                "ref": str,
                "alt": str,
                "idx": "int64",
                "allele_count": "int64",
            },
            keep_default_na=False,
            chunksize=chunk_rows,
        )
        for chunk in chunks:
            invalid = chunk["allele_count"].ne(2) | chunk[
                ["contig", "ref", "alt"]
            ].eq("").any(axis=1)
            if invalid.any():
                examples = chunk.loc[
                    invalid, ["contig", "position", "allele_count"]
                ].head(5)
                raise ValueError(
                    "Non-biallelic or incomplete variant rows found; refusing to build "
                    f"index. Examples: {examples.to_dict(orient='records')}"
                )

            table = pa.Table.from_pandas(
                chunk[ESSENTIAL_COLUMNS],
                preserve_index=False,
            )
            if writer is None:
                writer = pq.ParquetWriter(temporary_output, table.schema, compression="zstd")
            writer.write_table(table)
            row_count += table.num_rows

        if writer is None:
            raise ValueError("The exported variant table contains no rows.")
        finished_writer = writer
        writer = None
        finished_writer.close()
        os.replace(temporary_output, output_path)
        return row_count
    finally:
        if writer is not None:
            writer.close()
        if os.path.exists(temporary_output):
            os.remove(temporary_output)


def spark_conf() -> dict[str, str]:
    """Spark settings for the local-mode Hail session. The driver address is the loopback address: a cluster node whose
    hostname has an underscore (cn17_mgmt) otherwise gives Spark an invalid URL (spark://HeartbeatReceiver@cn17_mgmt:46155)
    and Hail fails to start. GWASFM_SPARK_DRIVER_HOST overrides it."""
    host = os.environ.get("GWASFM_SPARK_DRIVER_HOST", "127.0.0.1")
    return {
        "spark.driver.memory": os.environ.get("GWASFM_HAIL_DRIVER_MEMORY", "2g"),
        "spark.driver.host": host,
        "spark.driver.bindAddress": host,
        "spark.jars.packages": f"org.apache.hadoop:hadoop-aws:{HADOOP_AWS_VERSION}",
        "spark.hadoop.fs.s3a.aws.credentials.provider": "org.apache.hadoop.fs.s3a.AnonymousAWSCredentialsProvider",
        "spark.hadoop.fs.s3a.endpoint": "s3.us-east-1.amazonaws.com",
    }


def build_index(population: str, output_path: str, region: str | None = None) -> int:
    import hail as hl

    if population not in POPULATIONS:
        raise ValueError(f"Unsupported Pan-UKB population: {population}")

    output = Path(output_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing index: {output}")
    export_path = Path(f"{output_path}.export.tsv.bgz")
    temporary_output = Path(f"{output_path}.tmp")
    if export_path.exists() or temporary_output.exists():
        raise FileExistsError(f"Refusing to overwrite existing temporary index files for: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    tools_root = Path(
        os.environ.get(
            "GWASFM_TOOLS_ROOT",
            str(Path.home() / ".cache" / "gwasfinemapping"),
        )
    )
    log_dir = tools_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    table_path = (
        f"s3a://pan-ukb-us-east-1/ld_release/"
        f"UKBB.{population}.ldadj.variant.ht"
    )

    hl.init(
        master=os.environ.get("GWASFM_HAIL_MASTER", "local[2]"),
        log=str(log_dir / f"hail-panukb-{population}-{os.getpid()}.log"),
        spark_conf=spark_conf(),
    )
    try:
        table = hl.read_table(table_path)
        if "idx" not in table.row or "alleles" not in table.row:
            raise ValueError("Pan-UKB variant table must contain 'idx' and 'alleles' fields.")
        if region is not None:
            interval = hl.parse_locus_interval(region, reference_genome="GRCh37")
            table = hl.filter_intervals(table, [interval])
        table = table.annotate(
            contig=table.locus.contig,
            position=table.locus.position,
            ref=table.alleles[0],
            alt=table.alleles[1],
            allele_count=hl.len(table.alleles),
        )
        flat = table.key_by().select(*EXPORT_COLUMNS)
        flat.export(str(export_path))
        try:
            return write_index_from_export(str(export_path), output_path)
        finally:
            if export_path.exists():
                export_path.unlink()
    finally:
        hl.stop()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("population", choices=POPULATIONS)
    parser.add_argument("output", help="Output Parquet path")
    parser.add_argument("region", nargs="?", help="Optional GRCh37 region to limit the index")
    args = parser.parse_args()
    count = build_index(args.population, args.output, args.region)
    print(f"Wrote {count} variants to {args.output}")


if __name__ == "__main__":
    main()
