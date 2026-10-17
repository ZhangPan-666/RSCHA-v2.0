# Distribution build check, 2026-10-17

- Compared all 46 distributed C++ source/header files with the original OpenMP and MPI trees: SHA-256 identical.
- Built all six Windows executables using the distributed `scripts/build_windows.bat`: Work1, Work2, Work3 for OpenMP and MPI, all successful.
- Toolchain: local Intel oneAPI C++ 2025.3, Visual Studio C++ prerequisites, and Intel MPI.
- The table preparation script verified the canonical-residual table SHA-256 and copied it to both source trees before compilation.
- Compiler warnings include deprecated C runtime functions; there were no compilation failures.

This check establishes source completeness and local compilation. It does not replace numerical acceptance, exercise the six executables on scientific datasets, or certify another compiler/cluster environment. Linux build scripts require verification in the target environment.
