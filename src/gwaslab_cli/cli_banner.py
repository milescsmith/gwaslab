"""CLI mode banner on stderr; isolated so the launcher avoids heavy cli imports."""

from __future__ import annotations

import shlex
import sys
from importlib.metadata import PackageNotFoundError, version
from typing import List, Optional, TextIO

_BAR_LEN = 50


def _package_version() -> str:
    try:
        return version("gwaslab")
    except PackageNotFoundError:
        return "?.?.?"


def format_cli_mode_banner(prog: str, argv: List[str]) -> List[str]:
    """Return simple ASCII banner: rules, title, VERSION, COMMAND (one line), rules."""
    bar = "=" * _BAR_LEN
    ver = _package_version()
    title = "GWASLab CLI"
    cmd = shlex.join([prog] + list(argv))
    return [
        bar,
        title.center(_BAR_LEN),
        bar,
        f"VERSION: v{ver}",
        f"COMMAND: {cmd}",
        bar,
    ]


def emit_cli_mode_banner(
    prog: str,
    argv: List[str],
    *,
    quiet: bool,
    file: Optional[TextIO] = None,
) -> None:
    if quiet:
        return
    out = file if file is not None else sys.stderr
    print("\n".join(format_cli_mode_banner(prog, argv)) + "\n", file=out)
