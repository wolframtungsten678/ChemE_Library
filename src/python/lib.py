from __future__ import annotations

import ctypes
import os
from collections.abc import Iterator
from ctypes import c_double, c_int, c_int32
from ctypes.util import find_library
from pathlib import Path

_LIB_NAME = "ChemE_Library"
_ENV_OVERRIDE = "CHEME_LIBRARY_PATH"


def _library_candidates() -> Iterator[str | os.PathLike[str]]:
    """Yield potential locations for the compiled ChemE library."""
    env_override = os.environ.get(_ENV_OVERRIDE)
    if env_override:
        env_path = Path(env_override).expanduser()
        if env_path.is_file():
            yield env_path

    here = Path(__file__).resolve()
    repo_root = here.parents[2]
    build_dir = repo_root / "zig-out" / "lib"

    unix_name = f"lib{_LIB_NAME}"
    for suffix in (".so", ".dylib"):
        candidate = build_dir / f"{unix_name}{suffix}"
        if candidate.is_file():
            yield candidate

    windows_candidate = build_dir / f"{_LIB_NAME}.dll"
    if windows_candidate.is_file():
        yield windows_candidate

    located = find_library(_LIB_NAME)
    if located:
        yield located


def _load_library() -> ctypes.CDLL:
    errors: list[str] = []
    for candidate in _library_candidates():
        try:
            return ctypes.CDLL(os.fspath(candidate))
        except OSError as exc:
            errors.append(f"{candidate}: {exc}")

    joined = "\n  ".join(errors) if errors else "no candidates checked"
    raise FileNotFoundError(
        f"Unable to locate {_LIB_NAME} shared library. Tried:\n  {joined}"
    )


_LIB = _load_library()
_LIB.add.argtypes = (c_int32, c_int32)
_LIB.add.restype = c_int32


def add(a: int, b: int) -> int:
    """Return the sum computed by the Zig implementation."""
    return int(_LIB.add(c_int(a), c_int(b)))


class AddResult(ctypes.Structure):
    _fields_ = [("a", ctypes.c_int), ("b", ctypes.c_int), ("result", ctypes.c_int)]


_LIB.add2.argtypes = (c_int32, c_int32)
_LIB.add2.restype = AddResult


def add2(a: int, b: int) -> AddResult:
    """Return the sum computed by the Zig implementation."""
    result = _LIB.add2(c_int32(a), c_int32(b))
    return result


class SteamResult(ctypes.Structure):
    _fields_ = [
        ("ok", ctypes.c_bool),
        ("err_code", ctypes.c_int),
        ("phase_kind", ctypes.c_int),
        ("phase_region", ctypes.c_int),
        ("liquid_frac", ctypes.c_double),
        ("vapor_frac", ctypes.c_double),
        ("pressure", ctypes.c_double),
        ("temperature", ctypes.c_double),
        ("internal_energy", ctypes.c_double),
        ("enthalpy", ctypes.c_double),
        ("entropy", ctypes.c_double),
        ("cv", ctypes.c_double),
        ("cp", ctypes.c_double),
        ("speed_of_sound", ctypes.c_double),
        ("specific_volume", ctypes.c_double),
    ]


_LIB.getSteamEntryByPressureAndTemperature.argtypes = (c_double, c_double)
_LIB.getSteamEntryByPressureAndTemperature.restype = SteamResult


def getSteamEntryByPressureAndTemperature(
    pressure: float, temperature: float
) -> SteamResult:
    """Return the steam entry computed by the Zig implementation."""
    result = _LIB.getSteamEntryByPressureAndTemperature(
        c_double(pressure), c_double(temperature)
    )
    return result


__all__ = ["add", "add2", "getSteamEntryByPressureAndTemperature"]
