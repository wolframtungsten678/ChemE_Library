# tasks.py
import os
import pathlib
import sys
import sysconfig

import cffi


def _library_filename() -> str:
    if sys.platform == "darwin":
        return "libChemE_Library.dylib"
    if sys.platform == "win32":
        return "ChemE_Library.dll"
    return "libChemE_Library.so"


# Build the CFFI Python bindings
print("Building CFFI Module")
ffi = cffi.FFI()

project_root = pathlib.Path(__file__).resolve().parent
header_path = project_root / "src/core/root.h"
if not header_path.exists():
    raise FileNotFoundError(f"Cannot find header file at {header_path}")
ffi.cdef(header_path.read_text())

library_path = project_root / "zig-out" / "lib" / _library_filename()
if not library_path.exists():
    raise FileNotFoundError(
        f"Cannot find compiled library at {library_path}. "
        "Run `zig build` before building the Python bindings."
    )

ffi.set_source(
    "ChemE_Library._native",
    '#include "root.h"',
    include_dirs=[str(header_path.parent)],
    library_dirs=[str(library_path.parent)],
    libraries=["ChemE_Library"],
    extra_link_args=["-Wl,-rpath,."],
)

package_dir = project_root / "src/python/ChemE_Library"
package_dir.mkdir(parents=True, exist_ok=True)

# target_path = package_dir / f"*"
# ffi.compile(verbose=True, target=str(target_path))
ffi.compile()
