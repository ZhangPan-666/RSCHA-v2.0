# Regional application models

Final stored R-SCHA model coefficients used for the Cascadia and North American
Midcontinent examples in *Stable computation of revised spherical cap harmonics*.
These are the original application model files, preserved byte for byte.

| Region | File | Cap centre (longitude, geocentric latitude) | Coefficient entries |
|---|---|---|---:|
| Cascadia | [Results_RSCHAModel.dat](Cascadia/Results_RSCHAModel.dat) | -127 degrees, 45.5 degrees | 144,124 |
| Midcontinent | [Results_RSCHAModel.dat](Midcontinent/Results_RSCHAModel.dat) | -89 degrees, 45.8 degrees | 144,124 |

## File format

Each file is ASCII text with one number per row: nine header parameters followed
by 144,124 coefficient entries, for 144,133 rows in total. Do not reorder the
coefficients or discard stored zero-order sine entries.

| Header row | Meaning | Unit |
|---:|---|---|
| 1 | Cap-centre longitude | degrees east |
| 2 | Cap-centre geocentric latitude | degrees north |
| 3 | Cap half-angle | degrees |
| 4 | Internal Legendre truncation index | integer |
| 5 | External Legendre truncation index | integer |
| 6 | Mehler-family truncation limit | integer |
| 7 | Zero-eigenvalue-family order limit | integer |
| 8 | Outer shell radius | metres |
| 9 | Inner shell radius | metres |

Both models use a 7.5-degree cap, truncation settings (378, 7, 3, 3), and
outer/inner radii of 6,921,200 / 6,356,200 metres. The useful comparison cap in the
paper is 6.5 degrees. These native coefficients use the software's family-specific
normalizations and reference-radius convention; they are not global spherical
harmonic Gauss coefficients.

The coefficient blocks follow the public forward driver's order: internal
Legendre cosine and sine blocks (72,010 entries each), external Legendre cosine
and sine blocks (36 each), Mehler cosine and sine blocks (12 each), and
zero-eigenvalue cosine and sine blocks (4 each). Within a Legendre block the
indices run k = 0,...,K and m = 0,...,k. Mehler indices run p = 1,...,3 and
m = 0,...,3; zero-eigenvalue orders run m = 0,...,3.

## Forward evaluation

Use either released `Work2_RSCHAForward` driver with its required lookup table.
Copy the chosen model to the driver's working directory with the unchanged
basename `Results_RSCHAModel.dat`, and supply `Parameters_Coordinates.dat` in the
format documented with that driver. The model already contains the geometry
and truncation header. See the [OpenMP driver](../../src/v2.0-openmp/README.md)
and [MPI driver](../../src/v2.0-mpi/README.md) documentation at repository level.

The two files' SHA-256 digests and geometry metadata are recorded in
[models.json](models.json). Their hashes match the stored model identities used
in the application's figure-generation records. This upload does not rerun the
inversions or establish their unrecorded termination histories. Processed
observational input datasets and full inversion execution logs are not included.
Original source-data portals are listed in the paper. Existing repository terms
are described in LICENSE and LICENSING.md at repository level.
