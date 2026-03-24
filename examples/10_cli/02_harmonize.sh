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
# Reference files (defaults match examples/10_cli/07_pipeline.sh):
#   REF_SEQ / REF_VCF — gzipped FASTA + VCF built from INPUT below by
#   build_ref_from_sumstats.py (NEA=REF at POS, EA=ALT, AF=EAF) so harmonize/rsID match.
#   Override REF_SEQ, REF_VCF, and REF_INFER to use your own references; skip generation by
#   exporting them before running (remove or comment out the build line if both are set).
#   OUT_REF : directory for generated realistic_ref.* (default: ../../test/output)
#
# Requires:
#   pip install gwaslab
#   python3 (for build_ref_from_sumstats.py); bgzip/tabix optional (for VCF index)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="../../test/raw/realistic_sumstats.tsv"
mkdir -p output

OUT_REF="${OUT_REF:-../../test/output}"
mkdir -p "${SCRIPT_DIR}/${OUT_REF}"
REF_PREFIX="${OUT_REF}/realistic_ref"
python3 "${SCRIPT_DIR}/build_ref_from_sumstats.py" --input "${SCRIPT_DIR}/${INPUT}" --prefix "${SCRIPT_DIR}/${REF_PREFIX}"
REF_SEQ="${REF_SEQ:-${SCRIPT_DIR}/${REF_PREFIX}.fasta.gz}"
REF_VCF="${REF_VCF:-${SCRIPT_DIR}/${REF_PREFIX}.vcf.gz}"
REF_INFER="${REF_INFER:-${SCRIPT_DIR}/${REF_PREFIX}.vcf.gz}"


echo "=== [1] Harmonize with --ref-seq only (check_ref + flip; no rsID or --ref-infer) ==="
gwaslab \
    --input   "$INPUT" \
    --harmonize \
    --ref-seq "$REF_SEQ" \
    --output  "output/02_harmonized.tsv"

echo ""
echo "=== [2] Harmonize + assign rsID (--ref-rsid-vcf; default per-variant assignment) ==="
gwaslab \
    --input        "$INPUT" \
    --harmonize \
    --ref-seq      "$REF_SEQ" \
    --ref-rsid-vcf "$REF_VCF" \
    --output       "output/02_harmonized_with_rsid.tsv"

echo ""
echo "=== [3] Harmonize + strand inference (--ref-infer; defaults here use the demo VCF, not a population panel) ==="
gwaslab \
    --input      "$INPUT" \
    --harmonize \
    --ref-seq    "$REF_SEQ" \
    --ref-infer  "$REF_INFER" \
    --output     "output/02_harmonized_infer.tsv"

echo ""
echo "=== [4] Harmonize + strand inference with stricter MAF thresholds (sumstats vs ref VCF; requires --ref-infer) ==="
gwaslab \
    --input              "$INPUT" \
    --harmonize \
    --ref-seq            "$REF_SEQ" \
    --ref-infer          "$REF_INFER" \
    --maf-threshold      0.45 \
    --ref-maf-threshold  0.45 \
    --output             "output/02_harmonized_maf.tsv"

echo ""
echo "=== [5] Harmonize + rsID with --sweep-mode (chromosome sweep; faster on large data than [2]) ==="
gwaslab \
    --input        "$INPUT" \
    --harmonize \
    --ref-seq      "$REF_SEQ" \
    --ref-rsid-vcf "$REF_VCF" \
    --sweep-mode \
    --output       "output/02_harmonized_sweep.tsv"

echo ""
echo "=== [6] Harmonize with --threads (parallel work inside harmonization, e.g. normalize / ref steps) ==="
gwaslab \
    --input     "$INPUT" \
    --harmonize \
    --ref-seq   "$REF_SEQ" \
    --threads   4 \
    --output    "output/02_harmonized_parallel.tsv"

echo ""
echo "=== [7] QC with --remove, then harmonize (bad variants dropped before harmonization) ==="
gwaslab \
    --input     "$INPUT" \
    --qc \
    --remove \
    --harmonize \
    --ref-seq   "$REF_SEQ" \
    --output    "output/02_qc_then_harmonized.tsv"

echo ""
echo "Done. Output files written under output/"
