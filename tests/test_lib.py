from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = PROJECT_ROOT / "src" / "python"

if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

try:
    from lib import add, add2
except FileNotFoundError as exc:
    pytest.skip(
        f"Shared library missing ({exc}). Run `zig build` before executing the tests."
    )


def test_add_round_trip() -> None:
    assert add(21, 21) == 42


def test_add2_round_trip() -> None:
    assert add2(21, 21).result == 42
