#!/usr/bin/env bash
set -euo pipefail
bash "$(dirname "$0")/prepare_table.sh"
make -C "$(dirname "$0")/../src/v2.0-openmp"
make -C "$(dirname "$0")/../src/v2.0-mpi"
