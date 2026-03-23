#!/usr/bin/env bash
# =============================================================================
# 03_liftover.sh — Genome Build Liftover
# =============================================================================
# Lift GWAS summary statistics from one genome build to another.
# GWASLab automatically runs QC before liftover if --qc is not specified.
#
# IMPORTANT: the FROM build must match the actual build of the input file.
# --liftover <FROM> <TO> sets the input build automatically from <FROM>.
# If the build is unknown, use --infer-build to detect it from HapMap3 SNPs.
#
# Usage:
#   bash 03_liftover.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

INPUT="../0_sample_data/toy_data/dirty_sumstats.tsv"
mkdir -p output

echo "=== [1] Liftover hg19 → hg38 ==="
gwaslab \
    --input   "$INPUT" \
    --liftover 19 38 \
    --output  "output/03_lifted_hg38.tsv"

echo ""
echo "=== [2] Liftover hg38 → hg19 ==="
# Replace INPUT with an actual hg38 file
#gwaslab \
#    --input   "/path/to/hg38_sumstats.tsv" \
#    --liftover 38 19 \
#    --output  "output/03_lifted_hg19.tsv"

echo ""
echo "=== [3] Build unknown — infer build, then liftover to hg38 ==="
# --infer-build detects hg19 or hg38 from HapMap3 SNP coordinates.
# Combine with --liftover to let GWASLab decide the FROM build automatically.
gwaslab \
    --input       "$INPUT" \
    --qc \
    --infer-build \
    --liftover 19 38 \
    --output      "output/03_inferred_lifted_hg38.tsv"

echo ""
echo "=== [4] QC + liftover hg19 → hg38 in one step ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --remove \
    --liftover 19 38 \
    --output  "output/03_qc_lifted_hg38.tsv"

echo ""
echo "=== [5] Liftover and save as bgzipped + tabix-indexed VCF ==="
gwaslab \
    --input   "$INPUT" \
    --liftover 19 38 \
    --bgzip \
    --tabix \
    --to-fmt  vcf \
    --output  "output/03_lifted_hg38.vcf"

echo ""
echo "=== [6] Liftover hg19 → hg38 and output HapMap3 variants only ==="
gwaslab \
    --input    "$INPUT" \
    --qc \
    --liftover 19 38 \
    --hapmap3 \
    --output   "output/03_lifted_hg38_hapmap3.tsv"

echo ""
echo "Done. Output files written to output/"
