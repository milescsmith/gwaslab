#!/usr/bin/env bash
# =============================================================================
# 02_harmonize.sh — Harmonization
# =============================================================================
# Align GWAS summary statistics against a reference genome sequence.
# Harmonization fixes strand orientation, ref/alt allele assignment, and
# optionally assigns rsIDs from a reference VCF/HDF5.
#
# Usage:
#   bash 02_harmonize.sh
#
# Reference files (download separately):
#   REF_SEQ  : GRCh37/hg19 reference FASTA, e.g. human_g1k_v37.fasta
#   REF_VCF  : dbSNP VCF for rsID assignment
#   REF_INFER: 1000G VCF for strand inference
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

INPUT="../0_sample_data/toy_data/to_harmonize.tsv"
mkdir -p output

# ---- Paths to reference files -----------------------------------------------
# Update these to point at your local copies, or use:
#   gwaslab path hg19_fa    (to resolve configured paths)
REF_SEQ="${REF_SEQ:-/path/to/hg19.fa}"
REF_VCF="${REF_VCF:-/path/to/dbsnp_hg19.vcf.gz}"
REF_INFER="${REF_INFER:-/path/to/1000G_hg19.vcf.gz}"


echo "=== [1] Harmonize with reference FASTA (allele alignment only) ==="
gwaslab \
    --input   "$INPUT" \
    --harmonize \
    --ref-seq "$REF_SEQ" \
    --output  "output/02_harmonized.tsv"

echo ""
echo "=== [2] Harmonize + assign rsID from VCF ==="
gwaslab \
    --input        "$INPUT" \
    --harmonize \
    --ref-seq      "$REF_SEQ" \
    --ref-rsid-vcf "$REF_VCF" \
    --output       "output/02_harmonized_with_rsid.tsv"

echo ""
echo "=== [3] Harmonize + strand inference with 1000G panel ==="
gwaslab \
    --input      "$INPUT" \
    --harmonize \
    --ref-seq    "$REF_SEQ" \
    --ref-infer  "$REF_INFER" \
    --output     "output/02_harmonized_infer.tsv"

echo ""
echo "=== [4] Harmonize with custom MAF thresholds ==="
gwaslab \
    --input              "$INPUT" \
    --harmonize \
    --ref-seq            "$REF_SEQ" \
    --maf-threshold      0.45 \
    --ref-maf-threshold  0.45 \
    --output             "output/02_harmonized_maf.tsv"

echo ""
echo "=== [5] Harmonize + sweep mode (more aggressive ambiguous-SNP removal) ==="
gwaslab \
    --input       "$INPUT" \
    --harmonize \
    --ref-seq     "$REF_SEQ" \
    --sweep-mode \
    --output      "output/02_harmonized_sweep.tsv"

echo ""
echo "=== [6] Harmonize with parallel threads ==="
gwaslab \
    --input     "$INPUT" \
    --harmonize \
    --ref-seq   "$REF_SEQ" \
    --threads   4 \
    --output    "output/02_harmonized_parallel.tsv"

echo ""
echo "=== [7] Run QC then harmonize in one command ==="
gwaslab \
    --input     "$INPUT" \
    --qc \
    --remove \
    --harmonize \
    --ref-seq   "$REF_SEQ" \
    --output    "output/02_qc_then_harmonized.tsv"

echo ""
echo "Done. Output files written to output/"
