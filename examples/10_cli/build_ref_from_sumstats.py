#!/usr/bin/env python3
"""
Build minimal reference FASTA + VCF from a GWAS TSV so --harmonize and --ref-rsid-vcf
match the same CHR/POS/alleles as the sumstats.

Expects columns: CHR, POS, EA, NEA, EAF (and optionally SNPID). REF at each site is NEA[0];
ALT is EA. FASTA sequence at POS is the reference allele (NEA).

Usage:
  python build_ref_from_sumstats.py --input ../../test/raw/realistic_sumstats.tsv \\
      --prefix ../../test/output/realistic_ref
Writes: {prefix}.fasta.gz, {prefix}.vcf.gz, and {prefix}.vcf.gz.tbi when bgzip/tabix exist.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
from collections import defaultdict


def _read_variants(path: str) -> list[dict]:
    rows: list[dict] = []
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if not row.get("CHR") or row.get("CHR", "").strip() == "":
                continue
            try:
                chrom = int(float(str(row["CHR"]).replace("chr", "").replace("Chr", "").replace("CHR", "")))
            except (TypeError, ValueError):
                continue
            try:
                pos = int(float(str(row["POS"]).replace(",", "")))
            except (TypeError, ValueError):
                continue
            nea = str(row.get("NEA", "")).strip().upper()
            ea = str(row.get("EA", "")).strip().upper()
            if len(nea) != 1 or len(ea) != 1:
                continue
            eaf_s = str(row.get("EAF", "0.5")).strip()
            try:
                eaf = float(eaf_s)
            except ValueError:
                eaf = 0.5
            eaf = max(0.01, min(0.99, eaf))
            snpid = (row.get("SNPID") or "").strip() or f"rs{chrom}_{pos}"
            rows.append(
                {
                    "CHR": chrom,
                    "POS": pos,
                    "EA": ea,
                    "NEA": nea,
                    "EAF": eaf,
                    "SNPID": snpid,
                }
            )
    return rows


def _fill_n(seq: list[str]) -> None:
    bases = ("A", "T", "G", "C")
    for i in range(len(seq)):
        if seq[i] == "N":
            seq[i] = bases[i % 4]


def _write_fasta(variants: list[dict], out_path: str) -> dict[int, str]:
    """Return dict chrom -> sequence string."""
    by_chr: dict[int, list[tuple[int, str]]] = defaultdict(list)
    for v in variants:
        by_chr[v["CHR"]].append((v["POS"], v["NEA"]))

    sequences: dict[int, str] = {}
    for chrom in sorted(by_chr.keys()):
        max_pos = max(p for p, _ in by_chr[chrom])
        seq_len = max_pos + 1000
        seq = ["N"] * seq_len
        for pos, nea in by_chr[chrom]:
            if 0 < pos <= seq_len:
                seq[pos - 1] = nea[0]
        _fill_n(seq)
        sequences[chrom] = "".join(seq)

    with open(out_path, "w") as f:
        for chrom in sorted(sequences.keys()):
            f.write(f">{chrom}\n")
            s = sequences[chrom]
            for i in range(0, len(s), 60):
                f.write(s[i : i + 60] + "\n")
    return sequences


def _write_vcf(variants: list[dict], sequences: dict[int, str], out_path: str) -> None:
    variants = sorted(variants, key=lambda x: (x["CHR"], x["POS"]))
    with open(out_path, "w") as f:
        f.write("##fileformat=VCFv4.3\n")
        for chrom in sorted(sequences.keys()):
            f.write(f"##contig=<ID={chrom},length={len(sequences[chrom])}>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency of ALT">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            chrom, pos = v["CHR"], v["POS"]
            rsid = v["SNPID"]
            if not rsid.startswith("rs"):
                rsid = f"rs{chrom}_{pos}"
            ref, alt = v["NEA"], v["EA"]
            af = v["EAF"]
            f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\tAF={af}\n")


def _gzip_plain_to(src: str, dst_gz: str) -> None:
    with open(src, "rb") as f_in:
        with gzip.open(dst_gz, "wb", compresslevel=6) as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(src)


def _compress_vcf_and_tabix(vcf_plain: str, vcf_gz: str) -> None:
    bgzip = shutil.which("bgzip")
    if bgzip:
        try:
            subprocess.run([bgzip, "-f", vcf_plain], check=True, capture_output=True)
            produced = vcf_plain + ".gz"
            if os.path.isfile(produced) and os.path.abspath(produced) != os.path.abspath(vcf_gz):
                shutil.move(produced, vcf_gz)
        except subprocess.CalledProcessError:
            _gzip_plain_to(vcf_plain, vcf_gz)
    else:
        _gzip_plain_to(vcf_plain, vcf_gz)

    tabix = shutil.which("tabix")
    if tabix and os.path.isfile(vcf_gz):
        try:
            subprocess.run([tabix, "-p", "vcf", vcf_gz], check=True, capture_output=True)
        except subprocess.CalledProcessError:
            pass


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Sumstats TSV (CHR, POS, EA, NEA, EAF)")
    ap.add_argument(
        "--prefix",
        required=True,
        help="Output path prefix (writes prefix.fasta.gz, prefix.vcf.gz)",
    )
    args = ap.parse_args()
    variants = _read_variants(args.input)
    if not variants:
        print("No valid SNP rows found; need CHR, POS, EA, NEA (single bases).", file=sys.stderr)
        return 1

    prefix = args.prefix
    out_dir = os.path.dirname(os.path.abspath(prefix)) or "."
    os.makedirs(out_dir, exist_ok=True)

    fa_plain = prefix + ".fasta"
    vcf_plain = prefix + ".vcf"
    fa_gz = prefix + ".fasta.gz"
    vcf_gz = prefix + ".vcf.gz"

    sequences = _write_fasta(variants, fa_plain)
    _write_vcf(variants, sequences, vcf_plain)
    _gzip_plain_to(fa_plain, fa_gz)

    _compress_vcf_and_tabix(vcf_plain, vcf_gz)

    print(f"Wrote {fa_gz}")
    print(f"Wrote {vcf_gz}")
    if os.path.isfile(vcf_gz + ".tbi"):
        print(f"Wrote {vcf_gz}.tbi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
