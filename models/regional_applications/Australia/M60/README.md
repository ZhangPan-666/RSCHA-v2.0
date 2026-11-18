# Australian lithospheric magnetic model for tensor synthesis

[Results_RSCHAModel.dat](Results_RSCHAModel.dat) contains the coefficients used in
*Analytical computation of the lithospheric magnetic gradient tensor from high-resolution revised spherical cap harmonic models*. The model jointly fits airborne total-field magnetic anomalies and
a satellite-derived field synthesized at spherical harmonic degrees 16–120.

The nine native header rows encode a cap centre of 133° E, 27° S (geocentric),
a modelling half-angle of 22°, truncations
`(Kint, Kext, P, M) = (610, 15, 3, 60)`, and outer and inner radii of
6,921,200 and 6,356,200 m. The header is followed by 374,692 real coefficients.
Both side branches retain azimuthal orders 0–60. The target half-angle is 21°.
Heights in the paper are relative to the 6,371,200 m reference sphere.

Use the repository's [OpenMP forward driver](../../../../src/v2.0-openmp/README.md)
or [MPI forward driver](../../../../src/v2.0-mpi/README.md). Copy the file to
the driver's working directory as `Results_RSCHAModel.dat`. Supply geocentric
longitude and latitude in degrees and radius in metres in
`Parameters_Coordinates.dat`. For comparison with the paper, request the
geographic north, east and down output.

The model file has SHA-256 `8c1b5a3bf39c693dc35681d3d6d043e8e86592f2598826fcc4e45d34ada4ae21`.
