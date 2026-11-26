# Australian lithospheric magnetic model

*Calculating magnetic gradient tensor fields of the lithosphere from airborne total magnetic intensity and satellite magnetic vector data by revised spherical cap harmonic analysis*.

The near-surface total-field grid uses degrees 16–90 of the satellite-derived spherical harmonic model for long-wavelength restoration. The same model supplies vector data at degrees 16–120 on the 300, 350, 400, 450 and 500 km layers. The regional solution uses a coefficient-norm penalty of 0.001. At iteration 1,200 its relative normal-equation residual is 0.000545. The RMSE over 2,213,244 valid near-surface grid inputs is 12.34 nT.

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

The [raw coefficient file](https://raw.githubusercontent.com/ZhangPan-666/RSCHA-v2.0/main/models/regional_applications/Australia/Results_RSCHAModel.dat) has SHA-256 `73e7bdf8b7abe8ffa98ab3700e42f0ee1098599ac1b2f9507df1de424dfa6bc7`.

## Source data

The near-surface compilation uses the [2019 Geoscience Australia TMI grid](https://doi.org/10.26186/5e9cf3f2c0f1d). Original satellite observations are available from [GFZ for CHAMP](https://isdc.gfz.de/champ-isdc/access-to-the-champ-data), [ESA for Swarm](https://swarm-diss.eo.esa.int/), and the [Macao Science Satellite-1 data portal](https://mss.must.edu.mo/data_download.html).

## Funding for the Australian study

This work was supported by a Major Project of the Gansu Joint Research Fund (Grant No. 26JRRA1245), the Talent Support Program of Lanzhou Resources & Environment Voc-Tech University (Grant No. X2026RC-37), the National Natural Science Foundation of China (Grant Nos. 42574077, 42250103, and 42674085), and the Natural Science Foundation of Hubei Province (Grant No. 2024AFE003).
