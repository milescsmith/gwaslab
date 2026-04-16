#!/usr/bin/env bash
# Run GWASLab CLI against simulated chrX/chrY sumstats + minimal FASTA/VCF.
# Echoes input paths, a preview of the input table, the CLI invocation, and output preview.
#
# Usage (from repo root):
#   ./test/sexchr/run_cli_sexchr_simulated.sh
#   bash test/sexchr/run_cli_sexchr_simulated.sh
#
# If `gwaslab` is not on PATH, uses `python3 -m gwaslab_cli.cli` and sets PYTHONPATH to repo `src/`.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"
SRC_ROOT="$REPO_ROOT/src"
BUNDLED="$HERE/data/xy_sumstats.tsv"

TMP="$(mktemp -d)"
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT

export PYTHONPATH="${SRC_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

# Build minimal_xy.fa / .vcf.gz / xy_sumstats.tsv under $TMP (simulate_ref.py is Python-only).
mapfile -t _paths < <(
  HERE="$HERE" TMP_BUILD="$TMP" python3 <<'PY'
import os, sys, importlib.util

here = os.environ["HERE"]
tmp = os.environ["TMP_BUILD"]
name = "gwaslab_sexchr_simulate_ref_cli"
spec = importlib.util.spec_from_file_location(name, os.path.join(here, "simulate_ref.py"))
mod = importlib.util.module_from_spec(spec)
sys.modules[name] = mod
spec.loader.exec_module(mod)
p = mod.build_xy_reference(tmp)
print(p.sumstats_tsv)
print(p.fasta_path)
print(p.vcf_gz)
PY
)
SUMSTATS="${_paths[0]}"
FASTA="${_paths[1]}"
VCF="${_paths[2]}"

echo "=== Sex-chr CLI smoke test (simulated XY data) ==="
echo
echo "Input (simulator output, same rows as bundled data/xy_sumstats.tsv):"
echo "  SUMSTATS: $SUMSTATS"
echo "  REF_FASTA: $FASTA"
echo "  REF_VCF:   $VCF"
echo "  BUNDLED:   $BUNDLED (reference only)"
echo
echo "--- Input sumstats (first 8 lines) ---"
head -n 8 "$SUMSTATS"
echo

OUT_TSV="$TMP/cli_sexchr_out.tsv"

if command -v gwaslab >/dev/null 2>&1; then
  CLI=(gwaslab)
else
  CLI=(python3 -m gwaslab_cli.cli)
fi

cmd=(
  "${CLI[@]}"
  --input "$SUMSTATS"
  --fmt auto
  --build 38
  --fix-chr-pos
  --basic-check
  --harmonize
  --ref-seq "$FASTA"
  --ref-rsid-vcf "$VCF"
  --threads 1
  --output "$OUT_TSV"
  --to-fmt gwaslab
  --tab-fmt tsv
  --no-gzip
)

echo "--- CLI command ---"
printf '%q ' "${cmd[@]}"
echo
echo
echo "--- CLI stdout/stderr ---"
set +e
"${cmd[@]}"
rc=$?
set -e
echo
echo "--- exit code: $rc ---"
echo

shopt -s nullglob
written=("$TMP"/cli_sexchr_out.tsv*.tsv)
shopt -u nullglob

if [[ ${#written[@]} -gt 0 ]]; then
  out_file="${written[0]}"
  echo "--- Output file: $out_file (preview columns) ---"
  python3 - "$out_file" <<'PY'
import sys
import pandas as pd

path = sys.argv[1]
df = pd.read_csv(path, sep="\t")
want = ("CHR", "POS", "EA", "NEA", "STATUS", "rsID")
cols = [c for c in want if c in df.columns]
if cols:
    print(df[cols].to_string(index=False))
else:
    print(df.head().to_string(index=False))
PY
else
  echo "(No .tsv output found under $TMP matching cli_sexchr_out.tsv*)"
fi

exit "$rc"
