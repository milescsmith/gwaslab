"""Shared fixtures for CLI tests."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

# Ensure `src/` is importable when tests run without an editable install.
REPO_ROOT = Path(__file__).resolve().parents[2]
_SRC = REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

EXAMPLES_CLI = REPO_ROOT / "examples" / "10_cli"


@pytest.fixture(scope="session")
def gwaslab_bin(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Directory containing a `gwaslab` executable that runs the in-repo CLI."""
    bindir = tmp_path_factory.mktemp("gwaslab_bin")
    wrapper = bindir / "gwaslab"
    wrapper.write_text(
        f'#!/usr/bin/env bash\nexec "{sys.executable}" -m gwaslab_cli.main "$@"\n',
        encoding="utf-8",
    )
    wrapper.chmod(0o755)
    return bindir


@pytest.fixture
def gwaslab_env(gwaslab_bin: Path) -> dict[str, str]:
    """Environment: `gwaslab` on PATH, headless matplotlib, in-repo packages on PYTHONPATH."""
    env = os.environ.copy()
    env["PATH"] = str(gwaslab_bin) + os.pathsep + env.get("PATH", "")
    env.setdefault("MPLBACKEND", "Agg")
    # Prefer `src/` over any site-packages `gwaslab` so CLI tests match this checkout.
    src = str(REPO_ROOT / "src")
    env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
    return env


@pytest.fixture
def example_cli_workdir(tmp_path: Path) -> Path:
    """
    Mirror repo layout so scripts in examples/10_cli resolve ../../test/...:
      <tmp>/examples/10_cli/*.sh
      <tmp>/test -> <repo>/test
    """
    work = tmp_path / "examples" / "10_cli"
    work.mkdir(parents=True)
    for sh in EXAMPLES_CLI.glob("*.sh"):
        shutil.copy2(sh, work / sh.name)
    # Some example scripts call local helper Python utilities.
    for py in EXAMPLES_CLI.glob("*.py"):
        shutil.copy2(py, work / py.name)
    test_link = tmp_path / "test"
    if not test_link.exists():
        test_link.symlink_to(REPO_ROOT / "test", target_is_directory=True)
    return work
