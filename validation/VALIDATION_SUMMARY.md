# Public release validation summary

The acceptance record below belongs to the frozen numerical snapshot. This distribution also contains later application and compiler-compatibility changes; the historical acceptance does not certify all current distribution files. Current source hashes are recorded separately in PACKAGE_SHA256SUMS.txt.

The frozen acceptance package reports:

- primary arbitrary-precision comparisons: 40,995 value and angular-derivative requests;
- maximum primary scaled error: 1.11e-8 against a 1e-7 target;
- real-degree root refinement target: 1e-10, with maximum sampled relative Newton correction 2.83e-14;
- complete serial/three-rank MPI workflow: 20 inversion iterations, 70,264 valid forward coordinates and 7,136 coefficients;
- largest serial/MPI relative L2 difference: 1.81e-12;
- largest difference scaled by max(1, |a|, |b|): 3.56e-10;
- independent adjudication of the retained Bremer dispute queue: 249 distinct tuples, maximum new-method error 4.71e-9.

The raw Bremer comparison contains unavailable and disputed reference cases. Those cases are retained as provenance and are not counted as agreement. The finite grids and sampled parameter ranges do not certify a uniform error bound over the continuous parameter domain.

The local OpenMP/MPI synthetic regression covers 2, 4 and 8 OpenMP threads and 1, 3 and 6 MPI ranks, with seven 13-column output comparisons byte-identical. Full cluster-scale Intel 2019/1800-rank execution remains an environment-dependent check and is not claimed as a local certification.
