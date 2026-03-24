#!/usr/bin/env bash
# =============================================================================
# 01_report.sh — Generate GWASLab QC report using tutorial_v4 test data
# =============================================================================
# Use the same dataset/reference files as test/test_tutorial_v4.py and generate
# a GWASLab report (with harmonization and regional LD via VCF).
#
# Usage:
#   bash 01_report.sh
#
# Requires:
#   pip install gwaslab
#   python3
#   test/ref data files used by test/test_tutorial_v4.py
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT="${INPUT:-../../test/ref/bbj_t2d_hm3_chr7_variants.txt.gz}"
OUTDIR="${OUTDIR:-output}"
REPORT_HTML="${REPORT_HTML:-${OUTDIR}/tutorial_v4_report.html}"
SAVE_SUMSTATS="${SAVE_SUMSTATS:-${OUTDIR}/tutorial_v4_processed}"
REF_SEQ="${REF_SEQ:-../../test/ref/chr7.fasta.gz}"
REF_RSID_VCF="${REF_RSID_VCF:-../../test/ref/b157_2564.vcf.gz}"
REGIONAL_VCF="${REGIONAL_VCF:-../../test/ref/1kg_eas_hg19.chr7_126253550_128253550.vcf.gz}"

mkdir -p "${OUTDIR}"

INPUT_ABS="${SCRIPT_DIR}/${INPUT}"
REPORT_HTML_ABS="${SCRIPT_DIR}/${REPORT_HTML}"
SAVE_SUMSTATS_ABS="${SCRIPT_DIR}/${SAVE_SUMSTATS}"
REF_SEQ_ABS="${SCRIPT_DIR}/${REF_SEQ}"
REF_RSID_VCF_ABS="${SCRIPT_DIR}/${REF_RSID_VCF}"
REGIONAL_VCF_ABS="${SCRIPT_DIR}/${REGIONAL_VCF}"

echo ""
echo "============================================================"
echo "Generate report (HTML) with tutorial_v4 test data"
echo "============================================================"
gwaslab report \
    --input         "${INPUT_ABS}" \
    --fmt           gwaslab \
    --output        "${REPORT_HTML_ABS}" \
    --build         19 \
    --filter-region 7 1 300000000 \
    --harmonize \
    --ref-seq       "${REF_SEQ_ABS}" \
    --ref-rsid-vcf  "${REF_RSID_VCF_ABS}" \
    --ref-infer     "${REGIONAL_VCF_ABS}" \
    --ref-alt-freq  AF \
    --vcf           "${REGIONAL_VCF_ABS}" \
    --sig-level     5e-8 \
    --windowsizekb  500 \
    --save-sumstats "${SAVE_SUMSTATS_ABS}" \
    --save-fmt      gwaslab \
    --save-tab-fmt  tsv

echo ""
echo "Done."
echo "Report: ${REPORT_HTML_ABS}"
echo "Input: ${INPUT_ABS}"
echo "References: ${REF_SEQ_ABS}, ${REF_RSID_VCF_ABS}, ${REGIONAL_VCF_ABS}"
