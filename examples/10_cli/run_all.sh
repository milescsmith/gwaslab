#!/usr/bin/env bash
# =============================================================================
# run_all.sh — Run all example CLI scripts in order
# =============================================================================
# Usage (from anywhere):
#   bash examples/10_cli/run_all.sh
#   cd examples/10_cli && bash run_all.sh
#
# Skip 08_utility.sh (formatbook update + GWAS Catalog download need network):
#   SKIP_08=1 bash run_all.sh
#
# Requires: pip install gwaslab, python3 for build_ref_from_sumstats.py (02, 07).
#   09_variant_filters.sh uses the repo src CLI when ../../src/gwaslab_cli exists.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

SCRIPTS=(
  01_basic_qc.sh
  02_harmonize.sh
  03_liftover.sh
  04_format_conversion.sh
  05_plot.sh
  06_extract.sh
  07_pipeline.sh
  09_variant_filters.sh
)

if [[ "${SKIP_08:-}" != "1" ]]; then
  SCRIPTS+=(08_utility.sh)
else
  echo "SKIP_08=1 — skipping 08_utility.sh (network-heavy)"
fi

for s in "${SCRIPTS[@]}"; do
  echo ""
  echo "######################################################################"
  echo "# Running ${s}"
  echo "######################################################################"
  bash "./${s}"
done

echo ""
echo "######################################################################"
echo "# All selected scripts finished OK"
echo "######################################################################"
