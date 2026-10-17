# R-SCHA v2.0

Revised Spherical Cap Harmonic Analysis (R-SCHA), with OpenMP and MPI application drivers.

Repository: https://github.com/ZhangPan-666/RSCHA-v2.0

The release ZIP contains both complete source trees, the mathematical foundations handbook, validation records, build wrappers, and the canonical-residual lookup table. Download the ZIP from the repository's Releases page; GitHub's automatically generated source archives contain only files committed to the repository.

## Contents

- `src/v2.0-openmp/`: OpenMP drivers and headers, including the `SHRSCHA/` localization implementation.
- `src/v2.0-mpi/`: MPI drivers and headers, including `SHRSCHA/`.
- `docs/mathematics/`: mathematical foundations handbook PDF, corrections, and validation notes.
- `validation/`: historical acceptance record, validation scope, and checksums.
- `scripts/`: table preparation and build wrappers.
- `large_assets/cr.bin`: canonical-residual lookup table (158,806,776 bytes).

The applications are `Work1_RSCHAInversion` (regional inversion), `Work2_RSCHAForward` (field and tensor evaluation), and `Work3_LocalizeSHModel` (SH to R-SCHA localization). See each source directory's README for input files and cluster launch instructions. Measurement and global-model datasets are supplied separately.

## Build and run

Keep the extracted directory structure intact. The table preparation scripts verify `large_assets/cr.bin` against its SHA-256 and copy it into each source directory, as expected by the Makefiles and cluster scripts.

On Linux, load a C++14 compiler and a matching MPI C++ wrapper, then run:

```bash
bash scripts/build_linux.sh
```

For one variant, run `bash scripts/prepare_table.sh`, then `make -C src/v2.0-openmp` or `make -C src/v2.0-mpi`. Compiler selection can be overridden with `CXX`. The numerical build flags disable fast-math transformations.

On Windows, install Intel oneAPI C++ with its Visual Studio C++ prerequisites, then run `scripts\build_windows.bat`. Set `ONEAPI_ROOT` if oneAPI is installed outside its default directory. MPI builds also require Intel MPI (`I_MPI_ROOT`). The wrapper prepares the lookup table and builds all three applications for each available variant.

Run applications from their source directory or set `RSCHA_CR_TABLE_DIR` to an absolute directory containing `cr.bin`. Input/output directory arguments are described in the source READMEs. For Slurm, adjust partition, resource counts, compiler initialization, `WORK`, and `DATA_DIR` in `run.pbs` to match the target cluster before submission.

## Validation scope

`validation/acceptance.json` is a historical acceptance record for the frozen numerical implementation. It includes identity, arbitrary-precision comparison, roots and normalization, performance, independent Bremer-dispute adjudication, serial/MPI workflow, and projection checks. Those finite samples do not establish a uniform error theorem over the full real parameter domain.

This distribution also includes later application, progress-reporting, and compiler-compatibility changes. The historical acceptance record is not a new certification of every file in this distribution. The documented Work2 synthetic regression covers local OpenMP/MPI output consistency; full cluster-scale validation remains environment dependent. See `validation/VALIDATION_SUMMARY.md` and the source READMEs for scope.

## Licensing and external data

No software or document license has been selected for this release. Public availability alone does not grant a separate reuse license. Third-party model, geology, and coastline datasets are not included or relicensed here.
