"""Backward-compatible wrapper for moved CLI handlers."""

from gwaslab_cli.handlers import (
    load_sumstats,
    run_config,
    run_config_show,
    run_path,
    run_fb_list,
    run_fb_show,
    run_fb_update,
    run_version,
    run_download,
)

__all__ = [
    "load_sumstats",
    "run_config",
    "run_config_show",
    "run_path",
    "run_fb_list",
    "run_fb_show",
    "run_fb_update",
    "run_version",
    "run_download",
]
