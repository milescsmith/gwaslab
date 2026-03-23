#!/usr/bin/env bash
# =============================================================================
# 06_extract.sh — Extract Lead / Novel Variants
# =============================================================================
# Extract genome-wide significant loci or novel variants (not in GWAS Catalog).
#
# Usage:
#   bash 06_extract.sh
#
# Requires:
#   pip install gwaslab
#   Internet access for --extract novel (queries GWAS Catalog API)
# =============================================================================

set -euo pipefail

INPUT="../../test/raw/dirty_sumstats.tsv"
mkdir -p output

# --------------------------------------------------------------------------
# Lead variants
# --------------------------------------------------------------------------
echo "=== [1] Extract lead variants (default p < 5e-8, window 500 kb) ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --extract lead \
    --output  "output/06_lead_variants.tsv"

echo ""
echo "=== [2] Lead variants with relaxed p-value threshold (p < 1e-5) ==="
gwaslab \
    --input            "$INPUT" \
    --qc \
    --extract          lead \
    --sig-level-extract 1e-5 \
    --output           "output/06_lead_relaxed.tsv"

echo ""
echo "=== [3] Lead variants with narrower clumping window (250 kb) ==="
gwaslab \
    --input        "$INPUT" \
    --qc \
    --extract      lead \
    --windowsizekb 250 \
    --output       "output/06_lead_250kb.tsv"

# --------------------------------------------------------------------------
# Novel variants (requires GWAS Catalog API)
# --------------------------------------------------------------------------
echo ""
echo "=== [4] Extract novel variants for a trait (BMI, EFO_0004340) ==="
# The EFO ID specifies the trait used to pull known loci from GWAS Catalog
gwaslab \
    --input   "$INPUT" \
    --qc \
    --remove \
    --extract novel \
    --efo     EFO_0004340 \
    --output  "output/06_novel_bmi.tsv"

echo ""
echo "=== [5] Novel variants — return only truly novel hits ==="
gwaslab \
    --input      "$INPUT" \
    --qc \
    --remove \
    --extract    novel \
    --efo        EFO_0004340 \
    --only-novel \
    --output     "output/06_novel_only.tsv"

echo ""
echo "=== [6] Novel variants for multiple traits ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --remove \
    --extract novel \
    --efo     EFO_0004340 EFO_0004330 \
    --output  "output/06_novel_multi_trait.tsv"

echo ""
echo "Done. Output files written to output/"
