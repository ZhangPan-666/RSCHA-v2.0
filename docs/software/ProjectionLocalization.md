# Work3: Rotation and analytical projection from SH to R-SCHA

Follow the [build instructions](../../README.md#build-and-run) to prepare the lookup table and compile. Run the program from its source directory or set `RSCHA_CR_TABLE_DIR` to the absolute path of the directory containing `cr.bin`.

## Forward projection of two SH sources

```text
Work3_LocalizeSHModel.exe [input_directory] [output_directory]
```

Both directory arguments are optional. The input directory defaults to the current directory, and the output directory defaults to the input directory. Supply `Parameters_Modelling.dat` and `Parameters_SHModel.dat`. If `Parameters_SHModel_External.dat` is omitted, the external source is zero. Choose a separate output directory for each run to keep its results together.

Supply a complete triangular SH table from `n = 0` to `N`, sorted by degree and order, with columns `n`, `m`, `g`, and `h`. Use real Schmidt semi-normalization without the Condon–Shortley phase and a reference radius of 6,371,200 m. Magnetic synthesis starts at `n = 1`. The `n = 0` row is a placeholder, and the `m = 0` sine coefficient makes no contribution to the real potential. The internal and external inputs may have different maximum degrees. Pad them with zeros to the same `N`.

The nine modeling parameters are the cap-center longitude, latitude, and half-angle in degrees, the internal and external Legendre cutoffs `k`, the radial Mehler cutoff `p`, the azimuthal cutoff `m` shared by the Mehler and constant families, and the outer and inner sphere radii in meters, in that order.

The optional `Parameters_SideProjection.dat` accepts two lines with one nonnegative integer each. They override rows 6 and 7 of the modeling parameters. Each Mehler cosine or sine block contains `P_Max*(M_Max+1)` entries, ordered first by radial index. Each constant-family block contains `M_Max+1` entries. `Projection_Metadata.txt` records the cutoffs used in the run.

`Results_RSCHAModel.dat` retains the nine-line header and the four-family `G`/`H` order read by Work2. The program also writes `Projection_Metadata.txt` and `Results_RSCHASpectrum.dat`. The spectrum covers the two Legendre families.

The implementation uses degree-by-degree passive Wigner rotation, four source cross-blocks (`ii`, `ie`, `ei`, and `ee`), and both side-boundary families. The Mehler and constant families include a factor of `R_E/a`, which sets their scale relative to the earlier implementation. For nearly coincident degrees, integration over the cap yields the limiting value. The internal-source API in `LocalizeSH2RSCHA.h` calls the same implementation. Supply a complete triangular table, set `n_Start=1`, and use a root cache for the same cap half-angle and cutoffs.
