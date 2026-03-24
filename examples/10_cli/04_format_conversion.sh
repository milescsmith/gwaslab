#!/usr/bin/env bash
# =============================================================================
# 04_format_conversion.sh — Format Conversion
# =============================================================================
# Convert GWAS summary statistics between supported formats.
# Use 'gwaslab formatbook list' to see all available formats.
#
# Usage:
#   bash 04_format_conversion.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

INPUT="../../test/raw/dirty_sumstats.tsv"
mkdir -p output

# --------------------------------------------------------------------------
# List / inspect formats first
# --------------------------------------------------------------------------
echo "=== Available formats ==="
gwaslab formatbook list

echo ""
echo "=== Show column mapping for LDSC format ==="
gwaslab formatbook show ldsc

echo ""
echo "=== Show column mapping for VCF format ==="
gwaslab formatbook show vcf

# --------------------------------------------------------------------------
# Conversion examples
# --------------------------------------------------------------------------
echo ""
echo "=== [1] Convert to LDSC format (plain text) ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --to-fmt  ldsc \
    --no-gzip \
    --output  "output/04_converted.ldsc"

echo ""
echo "=== [2] Convert to LDSC format (gzipped) ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --to-fmt  ldsc \
    --output  "output/04_converted.ldsc.gz"

echo ""
echo "=== [3] Convert to VCF format ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --to-fmt  vcf \
    --no-gzip \
    --output  "output/04_converted.vcf"

echo ""
echo "=== [4] Convert to VCF + bgzip + tabix index ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --remove \
    --to-fmt  vcf \
    --bgzip \
    --tabix \
    --output  "output/04_converted.vcf"

echo ""
echo "=== [5] Convert to SSF (GWAS-SSF) format ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --to-fmt  ssf \
    --output  "output/04_converted.ssf.tsv"

echo ""
echo "=== [6] Convert to default gwaslab format (plain) ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --to-fmt  gwaslab \
    --no-gzip \
    --output  "output/04_converted.gwaslab.tsv"

echo ""
echo "=== [7] Output HapMap3 variants only ==="
gwaslab \
    --input   "$INPUT" \
    --qc \
    --hapmap3 \
    --to-fmt  ldsc \
    --output  "output/04_hapmap3.ldsc.gz"

echo ""
echo "Done. Output files written under output/"
