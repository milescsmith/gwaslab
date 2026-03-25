"""
Simulated GRCh38-style reference assets for chr X and chr Y (small contigs).

Writes:
  - ``xy_sumstats.tsv`` copy logic (see ``data/xy_sumstats.tsv`` in-repo)
  - ``minimal_xy.fa`` — chrX / chrY sequences with REF bases at variant POS = NEA
  - ``minimal_xy.vcf.gz`` + ``.tbi`` — same sites with AF in INFO (pysam tabix)
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Tuple

# (CHR label in sumstats, 1-based POS, NEA/ref, EA/alt, EAF, rsID, CHROM in VCF)
XY_VARIANTS: List[Tuple[str, int, str, str, float, str, str]] = [
    ("X", 100, "A", "G", 0.3, "rs900100", "chrX"),
    ("X", 101, "C", "T", 0.2, "rs900101", "chrX"),
    ("Y", 200, "G", "A", 0.25, "rs900200", "chrY"),
    ("Y", 201, "T", "C", 0.15, "rs900201", "chrY"),
]

# Contig lengths (must cover max POS)
LEN_CHR_X = 500
LEN_CHR_Y = 500


@dataclass
class XYReferencePaths:
    out_dir: str
    sumstats_tsv: str
    fasta_path: str
    vcf_gz: str

    @property
    def vcf_tbi(self) -> str:
        return self.vcf_gz + ".tbi"


def _build_fasta_sequences() -> Tuple[str, str]:
    xseq = ["N"] * LEN_CHR_X
    yseq = ["N"] * LEN_CHR_Y
    for chrom, pos, nea, _ea, _eaf, _rs, vcf_chrom in XY_VARIANTS:
        i = pos - 1
        if vcf_chrom == "chrX":
            xseq[i] = nea
        else:
            yseq[i] = nea
    return "".join(xseq), "".join(yseq)


def write_xy_sumstats(path: str) -> None:
    """Write TSV matching ``data/xy_sumstats.tsv``."""
    lines = [
        "SNPID\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN",
    ]
    betas = [0.01, 0.02, 0.03, 0.04]
    for idx, (chrom, pos, nea, ea, eaf, rs, _vcf) in enumerate(XY_VARIANTS):
        lines.append(
            f"{rs}\t{chrom}\t{pos}\t{ea}\t{nea}\t{eaf}\t{betas[idx]}\t0.05\t0.5\t1000"
        )
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def write_xy_fasta(path: str) -> None:
    xseq, yseq = _build_fasta_sequences()
    with open(path, "w", encoding="utf-8") as f:
        f.write(">chrX\n")
        for i in range(0, len(xseq), 60):
            f.write(xseq[i : i + 60] + "\n")
        f.write(">chrY\n")
        for i in range(0, len(yseq), 60):
            f.write(yseq[i : i + 60] + "\n")


def write_xy_vcf_and_index(vcf_gz_path: str) -> None:
    import pysam

    plain = vcf_gz_path.replace(".vcf.gz", ".vcf")
    with open(plain, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID=chrX,length={LEN_CHR_X}>\n")
        f.write(f"##contig=<ID=chrY,length={LEN_CHR_Y}>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="ALT AF">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _c, pos, nea, ea, eaf, rs, vch in XY_VARIANTS:
            f.write(f"{vch}\t{pos}\t{rs}\t{nea}\t{ea}\t.\tPASS\tAF={eaf}\n")

    pysam.tabix_compress(plain, vcf_gz_path, force=True)
    os.remove(plain)
    pysam.tabix_index(vcf_gz_path, preset="vcf", force=True)


def build_xy_reference(out_dir: str) -> XYReferencePaths:
    """Create all assets under ``out_dir`` (directory is created if missing)."""
    os.makedirs(out_dir, exist_ok=True)
    sumstats = os.path.join(out_dir, "xy_sumstats.tsv")
    fasta = os.path.join(out_dir, "minimal_xy.fa")
    vcf_gz = os.path.join(out_dir, "minimal_xy.vcf.gz")

    write_xy_sumstats(sumstats)
    write_xy_fasta(fasta)
    write_xy_vcf_and_index(vcf_gz)

    return XYReferencePaths(
        out_dir=out_dir,
        sumstats_tsv=sumstats,
        fasta_path=fasta,
        vcf_gz=vcf_gz,
    )


def bundled_sumstats_path() -> str:
    """Path to committed ``data/xy_sumstats.tsv`` (same rows as simulator)."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "data", "xy_sumstats.tsv")
