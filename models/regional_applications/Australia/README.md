# Australian lithospheric magnetic model

The inversion combines a regional airborne TMI
compilation with a satellite-derived field at degrees 16 to 120. The stored
solution uses a coefficient-norm penalty of 0.001 and 1,200 iterations. It
reached the iteration limit with a relative normal-equation residual of
7.87e-4. Its fit and spatial stability are reported in the manuscript.

[Results_RSCHAModel.dat](Results_RSCHAModel.dat) uses the single-column
R-SCHA format described in the [regional model catalog](../README.md).
The nine header rows give longitude 133 degrees east, geocentric latitude
27 degrees south, a 22-degree cap half-angle, truncation indices (610, 15,
3, 3), and outer and inner radii of 6,921,200 and 6,356,200 m. They are
followed by 374,236 coefficients in native order. The paper evaluates the
field within 21 degrees of the cap center and defines heights relative to
the 6,371,200 m sphere.

For forward evaluation, place the model in the working directory of
[Work2_RSCHAForward](../../../src/v2.0-openmp/README.md). Supply
`Parameters_Coordinates.dat` with longitude in degrees east, geocentric
latitude in degrees north, and geocentric radius in metres. Use the
geographic north, east and down output for comparison with the paper.

SHA-256 for the complete model file 598b1f624762089ebfc0fee9dc501b11b442157cc533b58b3ad6bbfcf2784f9b.
