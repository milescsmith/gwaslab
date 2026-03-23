"""
GWASLab CLI - Handler functions.

All helper functions used by the CLI live here so that package import remains
lightweight until actual commands require loading ``gwaslab``.
"""

from __future__ import annotations

import sys
import os
import json
from typing import Optional


def _resolve_path_keyword(gl, keyword: str) -> Optional[str]:
    """Resolve a path keyword from built-in paths or downloaded records."""
    if keyword in gl.options.paths:
        return gl.options.paths[keyword]
    return gl.get_path(keyword, verbose=False)


def load_sumstats(path: str, fmt: str, nrows: Optional[int] = None, build: str = "19"):
    """Load summary statistics from file.

    Parameters
    ----------
    path : str
        Path to input file.
    fmt : str
        Input format.
    nrows : int, optional
        Number of rows to read (for testing).
    build : str, optional
        Genome build of the input file (default: "19"). Encoded into the
        STATUS column at load time; must match the actual build of the input
        data.

    Returns
    -------
    gl.Sumstats
        Loaded Sumstats object.
    """
    import gwaslab as gl

    kwargs: dict = dict(fmt=fmt, build=build)
    if nrows is not None:
        kwargs["nrows"] = nrows
    return gl.Sumstats(path, **kwargs)


def run_config(args) -> None:
    """Handle the ``config`` subcommand."""
    import json
    import gwaslab as gl

    if args.json:
        print(json.dumps({"paths": gl.options.paths}, indent=2))
    else:
        for key, value in gl.options.paths.items():
            print(f"{key}: {value}")


def run_config_show(args) -> None:
    """Handle the ``config show`` subcommand."""
    import gwaslab as gl

    path = _resolve_path_keyword(gl, args.keyword)
    if not path:
        print(f"Path not found: {args.keyword}", file=sys.stderr)
        sys.exit(1)

    # For built-in JSON-backed config keys, display file content directly.
    if args.keyword in {"config", "reference", "formatbook"}:
        if not os.path.exists(path):
            print(f"Config file not found: {path}", file=sys.stderr)
            sys.exit(1)
        try:
            with open(path, "r", encoding="utf-8") as f:
                payload = json.load(f)
            print(json.dumps(payload, indent=2))
            return
        except Exception as exc:
            print(f"Failed to load JSON from {path}: {exc}", file=sys.stderr)
            sys.exit(1)

    # For non-JSON keys (e.g. data_directory or downloaded keyword), print resolved path.
    print(path)


def run_path(args) -> None:
    """Handle the ``path`` subcommand."""
    import gwaslab as gl

    path = _resolve_path_keyword(gl, args.keyword)
    if path:
        print(path)
    else:
        print(f"Path not found: {args.keyword}", file=sys.stderr)
        sys.exit(1)


def run_fb_list(args) -> None:
    """Handle the ``formatbook list`` subcommand."""
    import json
    import gwaslab as gl

    formats = gl.list_formats()
    if args.json:
        print(json.dumps({"formats": formats}, indent=2))
    else:
        for fmt in formats:
            print(fmt)


def run_fb_show(args) -> None:
    """Handle the ``formatbook show`` subcommand."""
    import json
    import gwaslab as gl

    mapping = gl.check_format(args.format)
    print(json.dumps({args.format: mapping}, indent=2))


def run_fb_update(args) -> None:
    """Handle the ``formatbook update`` subcommand."""
    import gwaslab as gl

    gl.update_formatbook()


def run_version(args) -> None:
    """Handle the ``version`` subcommand."""
    import gwaslab as gl

    gl.show_version()


def run_download(args) -> None:
    """Handle the ``download-sumstats`` subcommand."""
    import gwaslab as gl

    gl.download_sumstats(
        gcst_id=args.gcst_id,
        output_dir=args.output_dir,
        verbose=True,
    )

