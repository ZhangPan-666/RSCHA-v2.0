# R-SCHA v2.0 MPI instructions

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

The script retains its historical `run.pbs` name. It requests the `cpu` partition with `--time=0`, subject to cluster and account limits. The default allocation is 40 nodes with 45 MPI ranks per node, one thread per rank, and 170 GiB per node. The script verifies `cr.bin` and the required inputs, builds the selected application, and writes the results to `DATA_DIR`.

## Application inputs and outputs

- **Work1** uses the initial model, weights, and prior-model files. The default switches are `IncludeTMI=true`, `IncludeVector=false`, `IncludeDiffTMI=true`, and `IncludeDiffVector=true`. Enable `IncludeVector` when all three single-point vector input files are available. The script checks required files, then performs a one-process `--validate-only` pass before constructing the design matrix or evaluating the basis functions.
- **Work2** reads `Results_RSCHAModel.dat` and `Parameters_Coordinates.dat` and writes `Results_ForwardResults.dat`. The [Australian model example](../../models/regional_applications/Australia/README.md) gives the coefficient file and coordinate convention.
- **Work3** reads `Parameters_Modelling.dat` and `Parameters_SHModel.dat`. An external SH model and a side-projection configuration can also be supplied. It writes `Results_RSCHAModel.dat`, the spectrum, and projection metadata. See the [projection guide](../../docs/software/ProjectionLocalization.md).

## Compiler and runtime setup

The submission script loads `/opt/intel/bin/compilervars.sh intel64` when available, probes the C++ compiler, and passes it to `make`. The MPI variant also loads `mpivars.sh` from the Intel MPI 2019.7.217 installation when available. The C++ compilers are `icpc` or `icpx` for OpenMP and `mpiicpc` or `mpiicpx` for Intel MPI. GNU toolchains are supported. Use an MPI wrapper and runtime from the same installation.

The submission log records the compiler path and version. MPI normally launches through `mpirun`. If the cluster requires `srun`, configure `MPI_LAUNCHER` and the matching PMI or PMIx runtime. The local `cr.bin` is sufficient at runtime. Compiler compatibility helpers cover older C++14 toolchains, including alternatives to `std::clamp` and `std::filesystem` and path handling on Windows and POSIX. `make clean` can run before the compiler environment is loaded.

## Basis truncation

In `Parameters_Modelling.dat`, row 6 (`P_Max`) specifies the radial Mehler cutoff. Row 7 (`M_Max`) specifies the azimuthal cutoff shared by the Mehler and constant families. The Mehler modes use independent ranges `p = 1, …, P_Max` and `m = 0, …, M_Max`. Constant modes use `m = 0, …, M_Max`. Work3 accepts two nonnegative integers in `Parameters_SideProjection.dat` to override these rows in that order. The [user manual](../../docs/software/rscha_user_manual.pdf) gives the coefficient-block lengths and file layout.

## Fortran-style MPI inversion architecture (2026-10-15)

The MPI sum helpers retain their public names but use MPI_Gatherv, explicit rank-order summation on rank 0, and MPI_Bcast. The numerical sum helpers do not call MPI_Allreduce or MPI_Reduce; Work2 uses an integer-only MPI_Allreduce for progress counters. Root sum scratch storage is bounded to 64 MiB using blocks; counts and displacements remain within the legacy MPI int range. The design matrix stays local to each rank. Final predictions use indexed MPI_Gatherv and MPI_Bcast, with validation of complete, unique global row ownership.

Work1 distributes each data family in contiguous blocks, including the three components of each differenced vector pair. Kernel completion messages follow all-rank barriers. Rank-0 local progress and elapsed-time stage messages are flushed immediately; each rank also emits a final kernel-ready marker.

For fixed data weights (RobustSwitch=0), the preconditioner and RHS are constructed once and BiCGSTAB carries its shadow residual, current residual, search direction, alpha, omega, and beta across iterations. A recursively converged residual is verified with the actual weighted normal operator. When robust weights change, the Krylov state is restarted and the weighted system is rebuilt. Sigma, Huber/Tukey choices, damping, prior model, inactive coefficients, input filenames, output layouts, and solver function signature are retained. Iteration trajectories can differ from the previous one-step restart solver; historical numerical acceptance records remain tied to their original source snapshots.

run.pbs still selects WORK=1/2/3 and DATA_DIR. Work1 now performs the existing one-rank input validation, followed by --mpi-check on the full allocation with a 240-second communication timeout, before any root calculation or design matrix construction. --mpi-check is an optional diagnostic argument requiring no input data or atlas. A failed or timed-out check aborts the submission script. WORK=2/3 retain their interfaces and use the new sum helpers.

Work2 gathers only each rank's owned entries for the ten VBT fields, instead of collecting sparse full-length vectors from every rank. Geographic coordinate filtering, basis evaluation, tensor rotation, and the 13-column output format remain unchanged.

## Work2 forward progress

Rank 0 reports global progress during `forward local VBT calculation`. All ranks call `MPI_Allreduce` for the integer completion counts at the same local work rounds, so different rank speeds do not change the reduction order. The update preserves coordinate ownership, numerical results, MPI gathering, and the 13-column output format.

```text
[Forward progress] [##############--------------------------] 35.00% ...
```

About 25 updates per stage report the count, percentage, global point rate, and estimated remaining time. The progress bar covers local kernel evaluation. The subsequent `forward owned-value gather and broadcast` step retains stage logging. A long-running individual evaluation can delay the next synchronized update.

The October 17, 2026 update consolidated basis initialization, filters coordinates in one in-place pass without changing order or thresholds, and buffers output. The ten `LocalCart` arrays are allocated only on rank 0. For 474,496 valid coordinates this saves about 36.2 MiB per non-root rank (80 bytes per coordinate).

## Validation scope

Local Intel oneAPI C++14/OpenMP regressions used 2, 4, and 8 OpenMP threads and 1, 3, and 6 MPI ranks. The synthetic cases covered 3, 41, and 96 valid coordinates, all three orderings of the internal and external cutoffs, absent Mehler modes, invalid-coordinate filtering, and empty rank assignments. The seven comparisons produced byte-identical 13-column files and reached 100% progress. These checks cover the small local runs. The Intel 2019 cluster environment and full 40-node, 1,800-rank allocation require a separate cluster run. Historical acceptance records document their original source snapshots.
