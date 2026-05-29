"""Pytest configuration for bathy2strat."""

from __future__ import annotations

from pathlib import Path
import sys


def pytest_configure() -> None:
    """Ensure the project code directory is importable in tests."""
    code_dir = Path(__file__).resolve().parents[1]
    if str(code_dir) not in sys.path:
        sys.path.insert(0, str(code_dir))
