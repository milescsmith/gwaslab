#!/usr/bin/env bash
# =============================================================================
# 05_plot.sh — Visualization
# =============================================================================
# Generate plots from GWAS summary statistics.
# Available plot types: manhattan, qq, mqq, regional, forest
#
# Usage:
#   bash 05_plot.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

INPUT="../0_sample_data/toy_data/dirty_sumstats.tsv"
mkdir -p output/plots

echo "=== [1] Manhattan plot ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   manhattan \
    --output "output/plots/05_manhattan.png"

echo ""
echo "=== [2] QQ plot ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   qq \
    --output "output/plots/05_qq.png"

echo ""
echo "=== [3] Manhattan + QQ combined (mqq) ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   mqq \
    --output "output/plots/05_mqq.png"

echo ""
echo "=== [4] Manhattan with custom significance level ==="
gwaslab \
    --input     "$INPUT" \
    --qc \
    --plot      manhattan \
    --sig-level 1e-5 \
    --output    "output/plots/05_manhattan_sig1e5.png"

echo ""
echo "=== [5] Manhattan with custom y-axis limits ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   manhattan \
    --ylim   0 50 \
    --output "output/plots/05_manhattan_ylim.png"

echo ""
echo "=== [6] Regional plot (chromosome 6, 26–34 Mb) ==="
# Adjust --chr / --start / --end to a region with signal in your data
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   regional \
    --chr    6 \
    --start  26000000 \
    --end    34000000 \
    --output "output/plots/05_regional_chr6.png"

echo ""
echo "=== [7] Forest plot ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   forest \
    --output "output/plots/05_forest.png"

echo ""
echo "=== [8] Manhattan plot — quiet mode (suppress log messages) ==="
gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   manhattan \
    --quiet \
    --output "output/plots/05_manhattan_quiet.png"

echo ""
echo "Done. Plots written to output/plots/"
