#!/usr/bin/env bash
# =============================================================================
# 09_variant_filters.sh — PLINK-style filters and optional fix_* steps
# =============================================================================
# Demonstrates:
#   --fix-chr-pos, --extract / --exclude (ID lists), --extract-bed / --exclude-bed,
#   --chr, --maf / --max-maf, --mac, --snps-only
#
# --min-info is not run here (demo sumstats has no INFO column); see docs/CLI.md.
#
# Usage:
#   bash 09_variant_filters.sh
#
# Requires:
#   pip install gwaslab (recent version with variant-filter flags), **or**
#   run from a GWASLab git checkout — this script prepends ../../src to PYTHONPATH
#   and uses `python3 -m gwaslab_cli.cli` when that tree is present.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
if [[ -d "$REPO_ROOT/src/gwaslab_cli" ]]; then
  export PYTHONPATH="$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
  GWASLAB=(python3 -m gwaslab_cli.cli)
else
  GWASLAB=(gwaslab)
fi

INPUT="../../test/raw/realistic_sumstats.tsv"
OUT="output/09_variant_filters"
mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Helper list files (one variant ID per line; same format as PLINK --extract)
# ---------------------------------------------------------------------------
KEEP_LIST="${OUT}/keep.snplist"
DROP_LIST="${OUT}/drop.snplist"
BED_KEEP="${OUT}/region_chr1.bed"

# Keep first two SNPIDs from the table (skip header)
tail -n +2 "$INPUT" | cut -f1 | head -n 2 > "$KEEP_LIST"

# Drop one variant that is *not* in the keep list (so --exclude is visible on full data)
echo "6:3265891_A_G" > "$DROP_LIST"

# UCSC BED: 0-based half-open — narrow interval around first chr1 variant (POS 10949462)
printf '%s\t%d\t%d\n' 1 10949460 10949464 > "$BED_KEEP"

echo "=== [1] fix_chr + fix_pos then write (--fix-chr-pos) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --fix-chr-pos \
  --output "${OUT}/01_fix_chr_pos.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [2] --extract keep.snplist (only listed SNPIDs) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --extract "$KEEP_LIST" \
  --output "${OUT}/02_extract_only.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [3] --exclude drop.snplist (remove one variant) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --exclude "$DROP_LIST" \
  --output "${OUT}/03_exclude_one.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [4] --chr 1 2 (autosome subset) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --chr 1 2 \
  --output "${OUT}/04_chr1_chr2.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [5] --maf 0.20 (drops lowest-EAF variant; uses EAF→MAF) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --maf 0.20 \
  --output "${OUT}/05_maf020.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [6] --mac 50000 (with N + EAF; MAC ≈ 2×N×MAF) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --mac 50000 \
  --output "${OUT}/06_mac50k.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [7] --snps-only (single-nucleotide EA/NEA) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --snps-only \
  --output "${OUT}/07_snps_only.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [8] --extract-bed (keep variants overlapping BED) ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --build 19 \
  --extract-bed "$BED_KEEP" \
  --output "${OUT}/08_extract_bed.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "=== [9] Combined: exclude list + MAF + chr filter ==="
"${GWASLAB[@]}" \
  --input "$INPUT" \
  --exclude "$DROP_LIST" \
  --maf 0.19 \
  --chr 1 2 3 6 7 \
  --output "${OUT}/09_combined.tsv" \
  --to-fmt gwaslab \
  --tab-fmt tsv \
  --no-gzip

echo ""
echo "Done. Outputs under ${OUT}/"
wc -l "${OUT}"/*.tsv | sed 's/^/  /'
