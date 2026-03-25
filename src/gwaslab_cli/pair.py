"""GWASLab pair CLI for two-sumstats workflows."""

from __future__ import annotations

import argparse
import sys
from types import SimpleNamespace
from typing import Optional


SUPPORTED_OPERATIONS = ("merge", "miami")


def _find_operation_token(argv: list[str]) -> Optional[tuple[str, int]]:
    """Find operation token position in argv."""
    for idx, tok in enumerate(argv):
        if tok in SUPPORTED_OPERATIONS:
            return tok, idx
    return None


def _build_common_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=True, prog="gwaslab pair")
    parser.add_argument("--input1", "-i1", required=True, help="Input sumstats path for trait 1")
    parser.add_argument("--input2", "-i2", required=True, help="Input sumstats path for trait 2")
    parser.add_argument("--fmt1", default="auto", help="Input format for trait 1 (default: auto)")
    parser.add_argument("--fmt2", default="auto", help="Input format for trait 2 (default: auto)")
    parser.add_argument("--build", default="19", help="Input genome build for both traits (default: 19)")
    parser.add_argument(
        "--sync-alleles",
        action="store_true",
        help="Run fix_allele() on both inputs before pairing",
    )
    parser.add_argument(
        "--keep-all-variants",
        action="store_true",
        default=True,
        help="Keep all variants in outer merge (default: enabled)",
    )
    parser.add_argument(
        "--inner-join",
        action="store_true",
        help="Keep only variants shared by both datasets",
    )
    parser.add_argument("--threads", "-t", type=int, default=1, help="Thread count (reserved)")
    parser.add_argument("--quiet", "-q", action="store_true", help="Less logging")
    parser.add_argument("--output", "-o", required=True, help="Output path")
    return parser


def _build_operation_parser(op: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, prog=f"gwaslab pair {op}")
    if op == "miami":
        parser.add_argument("--sig-level", type=float, default=5e-8, help="Significance level")
    elif op == "merge":
        parser.add_argument("--sep", default="\t", help="Output delimiter for merged table")
    return parser


def parse_pair_args(argv: list[str]) -> tuple[SimpleNamespace, str]:
    """Parse pair command arguments with operation token in flexible position."""
    found = _find_operation_token(argv)
    if found is None:
        raise SystemExit(
            "pair: missing <operation>. Supported operations: merge, miami.\n"
            "Example:\n"
            "  gwaslab pair --input1 a.tsv --input2 b.tsv --fmt1 auto --fmt2 auto "
            "--build 19 --sync-alleles miami --output out.png"
        )
    op, op_idx = found
    filtered = argv[:op_idx] + argv[op_idx + 1 :]

    common_parser = _build_common_parser()
    common_ns, remaining = common_parser.parse_known_args(filtered)

    op_parser = _build_operation_parser(op)
    op_ns, op_unknown = op_parser.parse_known_args(remaining)
    if op_unknown:
        raise SystemExit(f"pair {op}: unrecognized arguments: {' '.join(op_unknown)}")

    merged = vars(common_ns)
    merged.update(vars(op_ns))
    return SimpleNamespace(**merged), op


def run_pair(argv: Optional[list[str]] = None) -> None:
    """Run pair workflow with operation-specific behavior."""
    if argv is None:
        argv = sys.argv[1:]
    args, operation = parse_pair_args(argv)

    import gwaslab as gl

    s1 = gl.Sumstats(args.input1, fmt=args.fmt1, build=args.build)
    s2 = gl.Sumstats(args.input2, fmt=args.fmt2, build=args.build)
    verbose = not args.quiet

    # Ensure basic coordinate consistency for robust pairing.
    s1.fix_chr(verbose=verbose)
    s1.fix_pos(verbose=verbose)
    s2.fix_chr(verbose=verbose)
    s2.fix_pos(verbose=verbose)

    if args.sync_alleles:
        s1.fix_allele(verbose=verbose)
        s2.fix_allele(verbose=verbose)

    keep_all_variants = False if args.inner_join else bool(args.keep_all_variants)
    pair = gl.SumstatsPair(s1, s2, keep_all_variants=keep_all_variants, verbose=verbose)

    if operation == "merge":
        pair.data.to_csv(args.output, sep=args.sep, index=False)
        return

    if operation == "miami":
        pair.plot_miami(save=args.output, sig_level=args.sig_level, verbose=verbose)
        return

    raise SystemExit(f"pair: unsupported operation: {operation}")
