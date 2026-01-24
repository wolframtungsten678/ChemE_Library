from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = PROJECT_ROOT / "src" / "python"

if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

try:
    from lib import add, add2, getSteamEntryByPressureAndTemperature
except FileNotFoundError as exc:
    pytest.skip(
        f"Shared library missing ({exc}). Run `zig build` before executing the tests."
    )


def test_add_round_trip() -> None:
    assert add(21, 21) == 42


def test_add2_round_trip() -> None:
    assert add2(21, 21).result == 42


def test_getSteamEntryByPressureAndTemperature_bad_round_trip() -> None:
    result = getSteamEntryByPressureAndTemperature(21, 21)
    assert result.ok == False
    # to do convert err_code to enum
    assert result.err_code == 148


def test_getSteamEntryByPressureAndTemperature_good_round_trip() -> None:
    result = getSteamEntryByPressureAndTemperature(40e6, 473.15)
    assert result.ok == True
    assert result.pressure == 40e6
    assert result.temperature == 473.15
