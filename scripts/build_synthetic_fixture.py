#!/usr/bin/env python3
"""Build a deterministic offline reference + 2-GWAS synthetic fixture for test.sh and benchmark.sh."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import subprocess
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


def write_fasta_with_fai(fasta_path: Path, contigs: dict[str, int]) -> None:
    """Write deterministic repeating ACGT FASTA + .fai index."""
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    width = 80
    unit = b"ACGT" * (width // 4)
    fai_lines: list[str] = []
    offset = 0
    with fasta_path.open("wb") as out:
        for name, length in contigs.items():
            header = f">{name}\n".encode("ascii")
            out.write(header)
            offset += len(header)
            fai_lines.append(f"{name}\t{length}\t{offset}\t{width}\t{width + 1}")
            full_lines, rem = divmod(length, width)
            block = (unit + b"\n") * 1000
            q, r_lines = divmod(full_lines, 1000)
            for _ in range(q):
                out.write(block)
            if r_lines:
                out.write((unit + b"\n") * r_lines)
            offset += full_lines * (width + 1)
            if rem:
                tail = (unit[:rem]) + b"\n"
                out.write(tail)
                offset += len(tail)
    Path(f"{fasta_path}.fai").write_text("\n".join(fai_lines) + "\n", encoding="utf-8")


def base_at(pos_1based: int) -> str:
    return "ACGT"[(pos_1based - 1) % 4]


def alt_base(ref: str) -> str:
    return {"A": "G", "C": "T", "G": "A", "T": "C"}[ref]


def write_chain(chain_path: Path, contigs: dict[str, int], shift: int = 1000) -> None:
    """Write UCSC chain mapping GRCh38 pos -> GRCh37 (pos - shift) cleanly."""
    chain_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for idx, (chrom, length) in enumerate(contigs.items(), start=1):
        span = length - shift - 1000
        lines.append(
            f"chain 1000000 chr{chrom} {length} + {shift} {shift + span} "
            f"chr{chrom} {length} + 0 {span} {idx}"
        )
        lines.append(f"{span}")
        lines.append("")
    with gzip.open(chain_path, "wt", encoding="utf-8") as gz:
        gz.write("\n".join(lines) + "\n")


def ar1_corr(positions: list[int], rho_per_10kb: float = 0.82) -> np.ndarray:
    pos = np.asarray(positions, dtype=np.float64)
    dist = np.abs(pos[:, None] - pos[None, :]) / 10_000.0
    R = np.power(rho_per_10kb, dist)
    np.fill_diagonal(R, 1.0)
    return R


def z_to_lp(z: float) -> float:
    p = math.erfc(abs(z) / math.sqrt(2.0))
    if p <= 1e-300:
        return 300.0
    return -math.log10(p)


def write_bgzip_vcf(vcf_gz: Path, rows: list[dict]) -> None:
    vcf_gz.parent.mkdir(parents=True, exist_ok=True)
    raw_vcf = vcf_gz.with_suffix("")
    header = "\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=8000000>",
        "##contig=<ID=2,length=10000000>",
        "##contig=<ID=3,length=5000000>",
        "##contig=<ID=6,length=35500000>",
        '##INFO=<ID=HM,Number=1,Type=String,Description="Harmonisation status">',
        '##FORMAT=<ID=ES,Number=1,Type=Float,Description="Effect size">',
        '##FORMAT=<ID=SE,Number=1,Type=Float,Description="Standard error">',
        '##FORMAT=<ID=LP,Number=1,Type=Float,Description="-log10 P-value">',
        '##FORMAT=<ID=SS,Number=1,Type=Float,Description="Sample size">',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Effect allele frequency">',
        '##FORMAT=<ID=AFS,Number=1,Type=String,Description="AF source">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ])
    rows_sorted = sorted(rows, key=lambda r: (int(r["chrom"]), int(r["pos"])))
    with raw_vcf.open("w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        for r in rows_sorted:
            sample = f"{r['es']:.6f}:{r['se']:.6f}:{r['lp']:.6f}:{r['ss']}:{r['af']:.4f}:study"
            fh.write(
                f"{r['chrom']}\t{r['pos']}\t{r['id']}\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                f"HM=direct\tES:SE:LP:SS:AF:AFS\t{sample}\n"
            )
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(vcf_gz), str(raw_vcf)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(vcf_gz)], check=True)
    raw_vcf.unlink(missing_ok=True)


def build_fixture(root: Path) -> dict:
    ref_root = root / "reference"
    gwas_root = root / "gwas"
    shift = 1000
    contigs = {
        "1": 8_000_000,
        "2": 10_000_000,
        "3": 5_000_000,
        "6": 35_500_000,
    }
    write_fasta_with_fai(ref_root / "genome" / "GRCh37.fa", contigs)
    write_chain(ref_root / "genome" / "hg38ToHg19.over.chain.gz", contigs, shift=shift)

    # Gene annotations: P00001 (GENEA on chr1, active cis), P00002 (GENENULL on chr3, sub-threshold cis)
    annot_dir = ref_root / "annotation"
    annot_dir.mkdir(parents=True, exist_ok=True)
    (annot_dir / "gene_regions_GRCh38.tsv").write_text(
        "uniprot_id\tgene_symbol\tensembl_id\tchrom\tstart\tend\tgenome_build\tsource\tversion\tsource_rows\n"
        "P00001\tGENEA\tENSG000001\t1\t3000000\t3050000\tGRCh38\tsynthetic\t1\t1\n"
        "P00002\tGENENULL\tENSG000002\t3\t2000000\t2040000\tGRCh38\tsynthetic\t1\t1\n",
        encoding="utf-8",
    )

    # Construct variant blocks on GRCh38 (note: base at pos on GRCh38 == base at pos-1000 on GRCh37 since 1000 % 4 == 0!)
    blocks = {
        "chr1_cis": ("1", [2_980_000 + i * 4_000 for i in range(25)]),
        "chr2_merged": (
            "2",
            [3_500_000 + i * 5_000 for i in range(15)]
            + [5_300_000 + i * 5_000 for i in range(15)],  # 1.8 Mb away -> merges under ±1.5 Mb (3 Mb) rule!
        ),
        "chr3_subthresh": ("3", [2_010_000 + i * 5_000 for i in range(12)]),
        "chr6_mhc": ("6", [31_200_000 + i * 6_000 for i in range(20)]),  # inside core GRC & extended MHC
    }

    panel_contigs, panel_pos, panel_ref, panel_alt, panel_idx = [], [], [], [], []
    matrix_blocks = []
    block_slices: dict[str, tuple[int, int, np.ndarray, list[int]]] = {}
    cursor = 0
    for bname, (chrom, gwas_positions) in blocks.items():
        R_b = ar1_corr(gwas_positions)
        m = len(gwas_positions)
        block_slices[bname] = (cursor, cursor + m, R_b, gwas_positions)
        for p38 in gwas_positions:
            p37 = p38 - shift
            r_allele = base_at(p37)
            a_allele = alt_base(r_allele)
            panel_contigs.append(chrom)
            panel_pos.append(p37)
            panel_ref.append(r_allele)
            panel_alt.append(a_allele)
            panel_idx.append(cursor)
            cursor += 1
        matrix_blocks.append(R_b)

    total_m = cursor
    full_panel_R = np.zeros((total_m, total_m), dtype=np.float64)
    offset = 0
    for R_b in matrix_blocks:
        m = R_b.shape[0]
        full_panel_R[offset:offset + m, offset:offset + m] = R_b
        offset += m

    panukb_dir = ref_root / "panukb"
    panukb_dir.mkdir(parents=True, exist_ok=True)
    whole_pq = panukb_dir / "UKBB.EUR.variant_index.parquet"
    pq.write_table(
        pa.table({
            "contig": panel_contigs,
            "position": pa.array(panel_pos, type=pa.int64()),
            "ref": panel_ref,
            "alt": panel_alt,
            "idx": pa.array(panel_idx, type=pa.int64()),
        }),
        whole_pq,
    )
    Path(f"{whole_pq}.json").write_text(
        json.dumps({"population": "EUR", "region": "genome", "build": "GRCh37", "rows": total_m}),
        encoding="utf-8",
    )
    local_ld_npz = ref_root / "local_eur_panel.npz"
    np.savez(local_ld_npz, matrix=full_panel_R, idx=np.arange(total_m, dtype=np.int64))

    rng = np.random.default_rng(20261117)
    n_samples = 20_000
    se_val = 1.0 / math.sqrt(n_samples)

    def make_block_rows(bname: str, chrom: str, causal_z: dict[int, float], noise_sd: float = 0.25) -> list[dict]:
        _, _, R_b, gwas_positions = block_slices[bname]
        m = len(gwas_positions)
        beta_std = np.zeros(m)
        for c_idx, z_target in causal_z.items():
            beta_std[c_idx] = z_target
        expected_z = R_b @ beta_std + rng.normal(0.0, noise_sd, size=m)
        rows = []
        for j, p38 in enumerate(gwas_positions):
            r_allele = base_at(p38)
            a_allele = alt_base(r_allele)
            z_j = float(expected_z[j])
            rows.append({
                "chrom": chrom,
                "pos": p38,
                "id": f"rs_{chrom}_{p38}",
                "ref": r_allele,
                "alt": a_allele,
                "es": z_j * se_val,
                "se": se_val,
                "lp": z_to_lp(z_j),
                "ss": n_samples,
                "af": 0.25 + 0.005 * (j % 20),
            })
        return rows

    # GWAS 1: P00001_olink (protein GWAS: 2-signal cis on chr1 + 2-lead merged trans on chr2 + MHC hit on chr6)
    g1_rows = (
        make_block_rows("chr1_cis", "1", {3: 9.5, 18: 8.2}, noise_sd=0.15)
        + make_block_rows("chr2_merged", "2", {4: 8.8, 22: 7.9}, noise_sd=0.15)
        + make_block_rows("chr3_subthresh", "3", {5: 3.1}, noise_sd=0.15)
        + make_block_rows("chr6_mhc", "6", {8: 10.2}, noise_sd=0.15)
    )
    g1_dir = gwas_root / "sun_synthetic" / "P00001_olink"
    g1_vcf = g1_dir / "P00001_olink.vcf.gz"
    write_bgzip_vcf(g1_vcf, g1_rows)
    (g1_dir / "metadata.json").write_text(json.dumps({
        "study_id": "sun_synthetic",
        "output_id": "P00001_olink",
        "uniprot_id": "P00001",
        "phenotype": "GENEA plasma protein",
        "ancestry": "EUR",
        "coordinate_build": "GRCh38",
    }, indent=2) + "\n", encoding="utf-8")
    (g1_dir / "qc.json").write_text(json.dumps({
        "coordinate_build": "GRCh38",
        "status": "PASS",
    }, indent=2) + "\n", encoding="utf-8")

    # GWAS 2: colon_rectal_sim (non-protein trait sharing causal index 4 on chr2 for coloc.susie test!)
    g2_rows = (
        make_block_rows("chr1_cis", "1", {}, noise_sd=0.2)
        + make_block_rows("chr2_merged", "2", {4: 9.1}, noise_sd=0.15)
        + make_block_rows("chr3_subthresh", "3", {5: 2.9}, noise_sd=0.15)
    )
    g2_dir = gwas_root / "cancer_synthetic" / "colon_rectal_sim"
    g2_vcf = g2_dir / "colon_rectal_sim.vcf.gz"
    write_bgzip_vcf(g2_vcf, g2_rows)

    # Manifest TSV for batch mode testing
    manifest_tsv = root / "manifest.tsv"
    manifest_tsv.write_text(
        "study_id\tgwas_id\ttrait\tancestry\tgene_id\tvcf\n"
        f"sun_synthetic\tP00001_olink\tGENEA plasma protein\tEUR\tP00001\t{g1_vcf}\n"
        f"cancer_synthetic\tcolon_rectal_sim\tColon and rectal cancer simulated\tEUR\t\t{g2_vcf}\n",
        encoding="utf-8",
    )

    return {
        "reference_root": str(ref_root),
        "local_ld_panel": str(local_ld_npz),
        "gwas1_vcf": str(g1_vcf),
        "gwas1_metadata": str(g1_dir / "metadata.json"),
        "gwas1_qc": str(g1_dir / "qc.json"),
        "gwas2_vcf": str(g2_vcf),
        "manifest_tsv": str(manifest_tsv),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="Directory in which to build synthetic fixture")
    args = parser.parse_args()
    info = build_fixture(Path(args.out).resolve())
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
