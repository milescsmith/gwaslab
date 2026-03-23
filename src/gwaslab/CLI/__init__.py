"""Backward-compatible GWASLab CLI namespace.

The canonical CLI implementation now lives in ``gwaslab_cli``.
This package remains as a compatibility layer for existing imports.
"""

from gwaslab_cli.cli import main

__all__ = ["main"]
