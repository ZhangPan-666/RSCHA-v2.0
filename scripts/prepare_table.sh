#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
(cd "$root/large_assets" && sha256sum -c cr.bin.sha256)
for variant in v2.0-openmp v2.0-mpi; do
    cp -- "$root/large_assets/cr.bin" "$root/src/$variant/cr.bin"
done
printf '%s\n' 'Verified and prepared cr.bin for both source trees.'
