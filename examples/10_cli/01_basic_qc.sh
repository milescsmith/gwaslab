#!/usr/bin/env bash
# =============================================================================
# 01_basic_qc.sh — Basic Quality Control
# =============================================================================
# Run QC on a GWAS summary statistics file and save the cleaned output.
#
# Usage:
#   bash 01_basic_qc.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

INPUT="../0_sample_data/toy_data/dirty_sumstats.tsv"
OUTPUT="output/01_cleaned.tsv"

mkdir -p output

echo "=== [1] Minimal QC (status flags only, no removal) ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --output "$OUTPUT"

echo ""
echo "=== [2] QC + remove bad-quality variants ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --remove \
    --output "output/01_cleaned_removed.tsv"

echo ""
echo "=== [3] QC + remove duplicates ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --remove-dup \
    --output "output/01_cleaned_no_dup.tsv"

echo ""
echo "=== [4] QC already normalizes indels by default ==="
# --normalize is redundant when --qc is used; shown here for clarity
gwaslab \
    --input  "$INPUT" \
    --qc \
    --output "output/01_cleaned_normalized.tsv"

echo ""
echo "=== [5] Full QC pipeline (remove + remove-dup; normalize is implicit) ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --remove \
    --remove-dup \
    --output "output/01_cleaned_full.tsv"

echo ""
echo "=== [6] QC with gzip output ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --remove \
    --output "output/01_cleaned.tsv.gz"

echo ""
echo "=== [7] QC — read only first 1000 rows (quick test) ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --nrows  1000 \
    --output "output/01_cleaned_1000.tsv"

echo ""
echo "Done. Output files written to output/"
