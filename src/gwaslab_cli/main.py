"""Lightweight launcher and dispatch for GWASLab CLI."""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version
from typing import Optional


def _print_version() -> None:
    """Print installed package version without importing gwaslab."""
    try:
        print(version("gwaslab"))
    except PackageNotFoundError:
        print("gwaslab (version unavailable)")


def main(argv: Optional[list[str]] = None) -> None:
    """Entry point for the lightweight GWASLab launcher."""
    if argv is None:
        argv = sys.argv[1:]

    # Fast path for trivial commands.
    if len(argv) == 1 and argv[0] == "version":
        _print_version()
        return

    # Delegate to full CLI for original parser behavior and docs.
    from gwaslab_cli.cli import main as gwaslab_main
    gwaslab_main(argv)


if __name__ == "__main__":
    main()

