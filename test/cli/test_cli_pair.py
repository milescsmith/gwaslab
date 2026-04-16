"""Tests for the dedicated `gwaslab pair` parser."""

from __future__ import annotations

import pytest

from gwaslab_cli.pair import parse_pair_args


def test_parse_pair_args_operation_in_middle() -> None:
    args, op = parse_pair_args(
        [
            "--input1",
            "trait1.tsv",
            "--input2",
            "trait2.tsv",
            "--fmt1",
            "auto",
            "--fmt2",
            "auto",
            "--build",
            "19",
            "--sync-alleles",
            "miami",
            "--output",
            "out.png",
        ]
    )
    assert op == "miami"
    assert args.input1 == "trait1.tsv"
    assert args.input2 == "trait2.tsv"
    assert args.sync_alleles is True
    assert args.output == "out.png"


def test_parse_pair_args_merge_with_custom_sep() -> None:
    args, op = parse_pair_args(
        [
            "--input1",
            "a.tsv",
            "--input2",
            "b.tsv",
            "merge",
            "--sep",
            ",",
            "--output",
            "merged.csv",
        ]
    )
    assert op == "merge"
    assert args.sep == ","
    assert args.output == "merged.csv"


def test_parse_pair_args_requires_operation() -> None:
    with pytest.raises(SystemExit):
        parse_pair_args(
            [
                "--input1",
                "a.tsv",
                "--input2",
                "b.tsv",
                "--output",
                "out.tsv",
            ]
        )
