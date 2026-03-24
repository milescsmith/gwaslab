#!/usr/bin/env bash
# =============================================================================
# 03_liftover.sh — Genome Build Liftover
# =============================================================================
# Lift GWAS summary statistics from one genome build to another.
# If --qc is not specified, CLI runs fix_chr + fix_pos before liftover
# (full QC still requires explicitly adding --qc).
#
# IMPORTANT: the FROM build must match the actual coordinates in the file.
# --liftover <FROM> <TO> passes FROM/TO into liftover; load uses build metadata when needed.
# Use --infer-build to log/detect hg19 vs hg38 from HapMap3 SNPs, then set --liftover FROM TO
# to match that detection (the CLI does not ignore the FROM you pass to --liftover).
#
# Usage:
#   bash 03_liftover.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

# hg19-style coordinates (must match --liftover FROM when not using infer-build alone).
INPUT="../../test/raw/to_harmonize.tsv"
mkdir -p output

echo "=== [1] Liftover hg19 → hg38 ==="
gwaslab \
    --input   "$INPUT" \
    --liftover 19 38 \
    --output  "output/03_lifted_hg38.tsv"

echo ""
echo "=== [2] Liftover hg38 → hg19 (skipped — uncomment and set INPUT to a real hg38 sumstats file) ==="
# Replace INPUT with an actual hg38 file
#gwaslab \
#    --input   "/path/to/hg38_sumstats.tsv" \
#    --liftover 38 19 \
#    --output  "output/03_lifted_hg19.tsv"

echo ""
echo "=== [3] QC + infer-build + liftover 19→38 ==="
# --infer-build updates detected build in metadata; --liftover still uses FROM/TO from
# --liftover (here 19→38). Align the first liftover argument with your true coordinate build.
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
echo "Done. Output files written under output/"
