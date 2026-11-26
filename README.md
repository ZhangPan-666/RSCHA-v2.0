# R-SCHA v2.0

Revised Spherical Cap Harmonic Analysis (R-SCHA), with OpenMP and MPI application drivers.

Repository: https://github.com/ZhangPan-666/RSCHA-v2.0

Software: [MIT](LICENSE). Mathematical handbook and its accompanying notes:

[CC BY 4.0](LICENSES/CC-BY-4.0.txt). See [licensing scope](LICENSING.md).

Use the main branch for the current sources, [program user manual](docs/software/rscha_user_manual.pdf), and [mathematical handbook](docs/mathematics/README.md).
The [v2.0 release ZIP](https://github.com/ZhangPan-666/RSCHA-v2.0/releases/download/v2.0.0/GitHub_Release_v2.0_20261017.zip) contains the original distribution and its lookup table.

## Contents

- [OpenMP source and instructions](src/v2.0-openmp/README.md), including the `SHRSCHA/` localization implementation.
- [MPI source and instructions](src/v2.0-mpi/README.md), including `SHRSCHA/`.
- `docs/software/`: program user manual and Work3 projection guide.
- `docs/mathematics/`: mathematical handbook PDF, LaTeX source archive, and audit records.
- `validation/`: historical acceptance record, validation scope, and checksums.
- `scripts/`: table preparation and build wrappers.
- `large_assets/cr.bin.sha256`: checksum of the canonical-residual lookup table.

The applications are `Work1_RSCHAInversion` (regional inversion), `Work2_RSCHAForward` (field and tensor evaluation), and `Work3_LocalizeSHModel` (SH to R-SCHA localization). See each source directory's README for input files and cluster launch instructions. Measurement and global-model datasets are supplied separately.

## Build and run

Extract `large_assets/cr.bin` (158,806,776 bytes) from the v2.0 release ZIP and place it in the same relative path in the current checkout. The preparation scripts verify the table against its SHA-256 and copy it into each source directory.

On Linux, use a C++14 compiler and a matching MPI C++ wrapper.

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

Project-authored software and supporting material are licensed under [MIT](LICENSE). The mathematical handbook and its accompanying project-authored notes are licensed under [CC BY 4.0](LICENSES/CC-BY-4.0.txt). See [LICENSING.md](LICENSING.md) for scope and attribution. Third-party model, geology, and coastline datasets retain their own terms and are not relicensed here.

## Regional application models

The final [Cascadia and Midcontinent model coefficients](models/regional_applications/README.md)

used in the manuscript's applications are available in native forward-driver format,

with geometry metadata and SHA-256 digests.
