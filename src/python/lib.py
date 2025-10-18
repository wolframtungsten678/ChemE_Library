from __future__ import annotations

import ctypes
import os
from collections.abc import Iterator
from ctypes import c_int32
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
    return int(_LIB.add(c_int32(a), c_int32(b)))


__all__ = ["add"]
