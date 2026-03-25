#!/usr/bin/env bash
# =============================================================================
# 07_pipeline.sh — Full Analysis Pipeline
# =============================================================================
# Runs several pipeline examples (A–F) in one script: QC, harmonize, liftover, plot, extract.
# Each block is a separate gwaslab invocation sharing the same INPUT/refs.
#
# Usage:
#   bash 07_pipeline.sh
#
# Reference files (defaults):
#   REF_SEQ / REF_VCF — gzipped FASTA + VCF built from INPUT below by
#   build_ref_from_sumstats.py (NEA=REF at POS, EA=ALT, AF=EAF) so harmonize/rsID match.
#   Override REF_SEQ and REF_VCF to use your own references; skip generation by exporting
#   both before running (the script will still run unless you remove that line).
#   OUT_REF   : directory for generated realistic_ref.* (default: ../../test/output)
#
# Requires:
#   pip install gwaslab
#   python3 (for build_ref_from_sumstats.py); bgzip/tabix optional (for VCF index)
# =============================================================================

set -euo pipefail

# Realistic hg19-scale positions (liftover-friendly); use dirty_sumstats.tsv for QC torture tests
INPUT="../../test/raw/realistic_sumstats.tsv"
mkdir -p output/pipeline output/pipeline/plots

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_REF="${OUT_REF:-../../test/output}"
mkdir -p "${SCRIPT_DIR}/${OUT_REF}"
REF_PREFIX="${OUT_REF}/realistic_ref"
python3 "${SCRIPT_DIR}/build_ref_from_sumstats.py" --input "${SCRIPT_DIR}/${INPUT}" --prefix "${SCRIPT_DIR}/${REF_PREFIX}"
REF_SEQ="${REF_SEQ:-${SCRIPT_DIR}/${REF_PREFIX}.fasta.gz}"
REF_VCF="${REF_VCF:-${SCRIPT_DIR}/${REF_PREFIX}.vcf.gz}"


# =============================================================================
# Pipeline A: QC → Save cleaned sumstats
# =============================================================================
echo "============================================================"
echo "Pipeline A: QC and clean"
echo "============================================================"
gwaslab \
    --input   "$INPUT" \
    --qc \
    --remove \
    --remove-dup \
    --normalize \
    --output  "output/pipeline/A_cleaned.tsv.gz"

echo ""


# =============================================================================
# Pipeline B: QC → Harmonize → Save
# =============================================================================
echo "============================================================"
echo "Pipeline B: QC → Harmonize"
echo "============================================================"
gwaslab \
    --input     "$INPUT" \
    --qc \
    --remove \
    --harmonize \
    --ref-seq   "$REF_SEQ" \
    --output    "output/pipeline/B_harmonized.tsv.gz"

echo ""


# =============================================================================
# Pipeline C: QC → Liftover → LDSC format
# =============================================================================
echo "============================================================"
echo "Pipeline C: QC → Liftover hg19→hg38 → LDSC"
echo "============================================================"
gwaslab \
    --input   "$INPUT" \
    --qc \
    --remove \
    --liftover 19 38 \
    --to-fmt  ldsc \
    --output  "output/pipeline/C_hg38.ldsc.gz"

echo ""


# =============================================================================
# Pipeline D: QC → Harmonize → Liftover → VCF (bgzipped + tabix)
# =============================================================================
echo "============================================================"
echo "Pipeline D: QC → Harmonize → Liftover → VCF"
echo "============================================================"
gwaslab \
    --input     "$INPUT" \
    --qc \
    --remove \
    --harmonize \
    --ref-seq   "$REF_SEQ" \
    --liftover  19 38 \
    --to-fmt    vcf \
    --bgzip \
    --tabix \
    --output    "output/pipeline/D_hg38.vcf"

echo ""


# =============================================================================
# Pipeline E: QC → Plot manhattan → Extract lead
# (Two separate gwaslab calls sharing the same input)
# =============================================================================
echo "============================================================"
echo "Pipeline E: QC → Manhattan plot + Lead extraction"
echo "============================================================"

gwaslab \
    --input  "$INPUT" \
    --qc \
    --plot   manhattan \
    --output "output/pipeline/plots/E_manhattan.png"

gwaslab \
    --input   "$INPUT" \
    --qc \
    --get lead \
    --output  "output/pipeline/E_lead_variants.tsv"

echo ""


# =============================================================================
# Pipeline F: QC → Harmonize → rsID assignment → Save
# =============================================================================
echo "============================================================"
echo "Pipeline F: QC → Harmonize → Assign rsID"
echo "============================================================"
gwaslab \
    --input        "$INPUT" \
    --qc \
    --remove \
    --harmonize \
    --ref-seq      "$REF_SEQ" \
    --ref-rsid-vcf "$REF_VCF" \
    --output       "output/pipeline/F_with_rsid.tsv.gz"

echo ""
echo "All examples complete. Outputs under output/pipeline/"
