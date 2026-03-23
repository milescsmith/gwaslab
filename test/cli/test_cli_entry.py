"""Smoke tests for the `gwaslab` CLI entry point (Python API)."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from pathlib import Path

import pytest

from gwaslab_cli.main import main

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_main_version_exits_cleanly() -> None:
    """`gwaslab version` fast path should not raise."""
    main(["version"])


def test_cli_help_does_not_crash() -> None:
    """`--help` should print and exit via SystemExit(0)."""
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0


def test_cli_qc_writes_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Minimal QC run with output file (matches public CLI examples)."""
    raw = REPO_ROOT / "test" / "raw" / "dirty_sumstats.tsv"
    if not raw.is_file():
        pytest.skip(f"missing test data: {raw}")
    out = tmp_path / "cleaned.tsv"
    monkeypatch.chdir(tmp_path)
    main(
        [
            "--input",
            str(raw),
            "--qc",
            "--output",
            str(out),
            "--quiet",
        ]
    )
    assert out.is_file() or any(tmp_path.glob("cleaned*"))


def test_cli_formatbook_list() -> None:
    """Subcommand `formatbook list` should run."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        main(["formatbook", "list"])
    text = buf.getvalue()
    assert text.strip()
