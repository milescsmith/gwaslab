#!/usr/bin/env bash
# =============================================================================
# 08_utility.sh — Utility Subcommands
# =============================================================================
# Examples of GWASLab's utility subcommands:
#   version        — show installed version
#   config         — inspect configured paths
#   path           — resolve a named reference path
#   formatbook     — list / inspect / update format definitions
#   download-sumstats — download from GWAS Catalog by accession ID
#
# Usage:
#   bash 08_utility.sh
#
# Requires:
#   pip install gwaslab
# =============================================================================

set -euo pipefail

mkdir -p output

# --------------------------------------------------------------------------
# version
# --------------------------------------------------------------------------
echo "=== Show installed GWASLab version ==="
gwaslab version

echo ""

# --------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------
echo "=== Show configured paths (human-readable) ==="
gwaslab config

echo ""
echo "=== Show configured paths (JSON) ==="
gwaslab config --json

echo ""

# --------------------------------------------------------------------------
# path
# --------------------------------------------------------------------------
echo "=== Resolve path for 'hg19_fa' ==="
gwaslab path hg19_fa || true     # '|| true' prevents script exit if not set

echo ""
echo "=== Resolve path for 'hg38_fa' ==="
gwaslab path hg38_fa || true

echo ""

# --------------------------------------------------------------------------
# formatbook
# --------------------------------------------------------------------------
echo "=== List all known formats ==="
gwaslab formatbook list

echo ""
echo "=== List formats as JSON ==="
gwaslab formatbook list --json

echo ""
echo "=== Show column mapping for 'ldsc' format ==="
gwaslab formatbook show ldsc

echo ""
echo "=== Show column mapping for 'vcf' format ==="
gwaslab formatbook show vcf

echo ""
echo "=== Show column mapping for 'ssf' format ==="
gwaslab formatbook show ssf

echo ""
echo "=== Update formatbook from upstream ==="
gwaslab formatbook update

echo ""

# --------------------------------------------------------------------------
# download-sumstats
# --------------------------------------------------------------------------
echo "=== Download sumstats from GWAS Catalog (GCST90270926) ==="
# Files are downloaded to --output-dir (defaults to current directory)
gwaslab download-sumstats GCST90270926 --output-dir output/

echo ""
echo "Done."
