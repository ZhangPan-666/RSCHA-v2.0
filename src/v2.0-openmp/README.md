# R-SCHA v2.0 OpenMP instructions

This source directory contains the project headers and the `SHRSCHA/` modules. Use the [lookup-table preparation instructions](../../README.md#build-and-run) to place the verified `cr.bin` in this directory before deploying it to the cluster.

## Slurm submission

Edit the user settings at the top of `run.pbs`:

```bash
WORK=1                         # 1: inversion; 2: forward synthesis; 3: SH to R-SCHA
DATA_DIR="data/Work1"           # Input and output directory
```

The script runs one application at a time. For the second survey, set `DATA_DIR="data/Work2"`. If the data files are in the source directory, set `DATA_DIR="."`. Relative paths are resolved from the submission directory. Absolute paths and paths with spaces are supported.

From the directory containing the script, submit it with Slurm:

```bash
sbatch run.pbs
```

The script retains its historical `run.pbs` name. It requests the `cpu` partition with `--time=0`, subject to cluster and account limits. The default allocation is one node, one process, and 45 OpenMP threads. The script verifies `cr.bin` and the required inputs, builds the selected application, and writes the results to `DATA_DIR`.

## Application inputs and outputs

- **Work1** uses the initial model, weights, and prior-model files. The default switches are `IncludeTMI=true`, `IncludeVector=false`, `IncludeDiffTMI=true`, and `IncludeDiffVector=true`. Enable `IncludeVector` when all three single-point vector input files are available. The script checks required files, then performs a one-process `--validate-only` pass before constructing the design matrix or evaluating the basis functions.
- **Work2** reads `Results_RSCHAModel.dat` and `Parameters_Coordinates.dat` and writes `Results_ForwardResults.dat`. The [Australian model example](../../models/regional_applications/Australia/README.md) gives the coefficient file and coordinate convention.
- **Work3** reads `Parameters_Modelling.dat` and `Parameters_SHModel.dat`. An external SH model and a side-projection configuration can also be supplied. It writes `Results_RSCHAModel.dat`, the spectrum, and projection metadata. See the [projection guide](../../docs/software/ProjectionLocalization.md).

## Compiler and runtime setup

The submission script loads `/opt/intel/bin/compilervars.sh intel64` when available, probes the C++ compiler, and passes it to `make`. The MPI variant also loads `mpivars.sh` from the Intel MPI 2019.7.217 installation when available. The C++ compilers are `icpc` or `icpx` for OpenMP and `mpiicpc` or `mpiicpx` for Intel MPI. GNU toolchains are supported. Use an MPI wrapper and runtime from the same installation.

The submission log records the compiler path and version. MPI normally launches through `mpirun`. If the cluster requires `srun`, configure `MPI_LAUNCHER` and the matching PMI or PMIx runtime. The local `cr.bin` is sufficient at runtime. Compiler compatibility helpers cover older C++14 toolchains, including alternatives to `std::clamp` and `std::filesystem` and path handling on Windows and POSIX. `make clean` can run before the compiler environment is loaded.

## Basis truncation

In `Parameters_Modelling.dat`, row 6 (`P_Max`) specifies the radial Mehler cutoff. Row 7 (`M_Max`) specifies the azimuthal cutoff shared by the Mehler and constant families. The Mehler modes use independent ranges `p = 1, …, P_Max` and `m = 0, …, M_Max`. Constant modes use `m = 0, …, M_Max`. Work3 accepts two nonnegative integers in `Parameters_SideProjection.dat` to override these rows in that order. The [user manual](../../docs/software/rscha_user_manual.pdf) gives the coefficient-block lengths and file layout.

## Work2 forward progress

During the potential, field, and tensor kernel calculation, an atomic counter tracks completed coordinates. The OpenMP reporting block prints about 25 updates per stage without interleaving or reversing progress messages. The serial fallback uses the same coordinate count.

```text
[Forward progress] [##############--------------------------] 166074/474496 (35.00%); 120.0 points/s; ETA 42.8 min
```

Updates show the number of completed points, percentage, average rate, and estimated remaining time. Lines are flushed immediately to the Slurm log. Rotation and file output retain their own stage messages. A long-running individual evaluation can delay the next progress update.

Work2 uses dynamic scheduling and retains its 13-column format. The October 16, 2026 update consolidated basis initialization and progress callbacks, filters coordinates in one in-place pass without changing order or thresholds, and buffers output. To deploy an updated driver, keep the filename `Work2_RSCHAForward.cpp` and resubmit with `WORK=2`. The build picks up the source change for the new run.

## Validation scope

Local Intel oneAPI C++14/OpenMP regressions used 2, 4, and 8 OpenMP threads and 1, 3, and 6 MPI ranks. The synthetic cases covered 3, 41, and 96 valid coordinates, all three orderings of the internal and external cutoffs, absent Mehler modes, invalid-coordinate filtering, and empty rank assignments. The seven comparisons produced byte-identical 13-column files and reached 100% progress. These checks cover the small local runs. The Intel 2019 cluster environment and full 40-node, 1,800-rank allocation require a separate cluster run. Historical acceptance records document their original source snapshots.
