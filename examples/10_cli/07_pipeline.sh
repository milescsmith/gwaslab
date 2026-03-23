#!/usr/bin/env bash
# =============================================================================
# 07_pipeline.sh — Full Analysis Pipeline
# =============================================================================
# Chains QC → Harmonization → Liftover → Extract → Plot in a single script.
# Demonstrates real-world workflow patterns.
#
# Usage:
#   bash 07_pipeline.sh
#
# Reference files (update paths or set env variables):
#   REF_SEQ   : hg19 reference FASTA
#   REF_VCF   : dbSNP VCF for rsID annotation
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

INPUT="../0_sample_data/toy_data/dirty_sumstats.tsv"
mkdir -p output/pipeline output/pipeline/plots

REF_SEQ="${REF_SEQ:-/path/to/hg19.fa}"
REF_VCF="${REF_VCF:-/path/to/dbsnp_hg19.vcf.gz}"


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
    --extract lead \
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
echo "All pipelines complete. Output in output/pipeline/"
