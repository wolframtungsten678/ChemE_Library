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
    assert result.err_code == 0
    assert result.phase_kind == 2  # NonCritical
    assert result.phase_region == 0  # Liquid
    assert result.liquid_frac == 0
    assert result.vapor_frac == 0

    assert result.pressure == pytest.approx(40e6, abs=1e-6)
    assert result.temperature == pytest.approx(473.15, abs=1e-6)
    assert result.internal_energy == pytest.approx(825.228016170348e3, abs=1e-3)
    assert result.enthalpy == pytest.approx(870.124259682489e3, abs=1e-3)
    assert result.entropy == pytest.approx(2.275752861241e3, abs=1e-3)
    assert result.cv == pytest.approx(3.292858637199e3, abs=1e-3)
    assert result.cp == pytest.approx(4.315767590903e3, abs=1e-3)
    assert result.speed_of_sound == pytest.approx(1457.418351596083, abs=1e-3)
    assert result.specific_volume == pytest.approx(0.001122406088, abs=1e-9)
