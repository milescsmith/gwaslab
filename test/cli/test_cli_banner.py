"""CLI mode banner (simple stderr header + VERSION / COMMAND lines)."""

from __future__ import annotations

import io
import re
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(autouse=True)
def _ensure_src_path() -> None:
    src = str(REPO_ROOT / "src")
    if src not in sys.path:
        sys.path.insert(0, src)


def test_format_cli_mode_banner_layout() -> None:
    from gwaslab_cli.cli_banner import format_cli_mode_banner

    lines = format_cli_mode_banner("gwaslab", ["version", "--help"])
    bar = "=" * 50
    assert lines[0] == bar
    assert len(lines[1]) == 50
    assert lines[1].strip() == "GWASLab CLI"
    assert lines[2] == bar
    assert re.match(r"^VERSION: v[\d.?]+$", lines[3])
    assert lines[4] == "COMMAND: gwaslab version --help"
    assert lines[5] == bar


def test_format_cli_mode_banner_version_and_command_lines() -> None:
    from gwaslab_cli.cli_banner import format_cli_mode_banner

    lines = format_cli_mode_banner("gwaslab", ["format", "--input", "raw.txt.gz", "--output", "formatted.tsv.gz"])
    assert lines[3].startswith("VERSION: v")
    assert lines[4] == (
        "COMMAND: gwaslab format --input raw.txt.gz --output formatted.tsv.gz"
    )


def test_format_cli_mode_banner_command_single_line_long_args() -> None:
    from gwaslab_cli.cli_banner import format_cli_mode_banner

    argv = ["--input", "x" * 60, "--output", "out.tsv"]
    lines = format_cli_mode_banner("gwaslab", argv)
    assert lines[4].startswith("COMMAND: gwaslab ")
    assert "x" * 60 in lines[4]


def test_emit_cli_mode_banner_respects_quiet() -> None:
    from gwaslab_cli.cli_banner import emit_cli_mode_banner

    buf = io.StringIO()
    emit_cli_mode_banner("gwaslab", ["config"], quiet=True, file=buf)
    assert buf.getvalue() == ""


def test_emit_cli_mode_banner_writes_banner_to_file() -> None:
    from gwaslab_cli.cli_banner import emit_cli_mode_banner

    buf = io.StringIO()
    emit_cli_mode_banner("gwaslab", ["path", "config"], quiet=False, file=buf)
    text = buf.getvalue()
    assert "VERSION:" in text
    assert "COMMAND: gwaslab path config" in text
    assert text.count("=" * 50) >= 3


def test_cli_main_emits_cli_mode_banner_to_stderr() -> None:
    from gwaslab_cli.cli import main as cli_main

    err = io.StringIO()
    out = io.StringIO()
    with redirect_stderr(err):
        with redirect_stdout(out):
            cli_main(["version"])
    stderr = err.getvalue()
    assert "GWASLab CLI" in stderr
    assert "VERSION: v" in stderr
    assert "COMMAND: gwaslab version" in stderr
    assert "=" * 50 in stderr


def test_launcher_version_fast_path_emits_cli_mode_banner() -> None:
    from gwaslab_cli.main import main as launcher_main

    err = io.StringIO()
    out = io.StringIO()
    with redirect_stderr(err):
        with redirect_stdout(out):
            launcher_main(["version"])
    stderr = err.getvalue()
    assert "GWASLab CLI" in stderr
    assert "COMMAND: gwaslab version" in stderr


def test_cli_main_quiet_suppresses_banner() -> None:
    from gwaslab_cli.cli import main as cli_main

    err = io.StringIO()
    out = io.StringIO()
    with redirect_stderr(err):
        with redirect_stdout(out):
            cli_main(["--quiet", "config"])
    assert "COMMAND:" not in err.getvalue()
    assert "config" in out.getvalue().lower() or "paths" in out.getvalue().lower()
