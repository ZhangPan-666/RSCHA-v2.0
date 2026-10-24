# Manuscript reproducibility records - revision manuscript-20261024-r2

This supplement accompanies the revised manuscript. It is distributed with the
source package under GitHub docs/manuscript; the software release remains v1.0.
It contains mathematical request coordinates,
aggregate statistics, model settings and plotted spectral series; it contains no
raw magnetic observations or restricted geological layers.

## Numerical references

`numerical_validation/README.md` documents all 40,995 archived requests and the
portable Python reference-only entry point. A six-point smoke regeneration and
24 independent Python controls passed. The compact historical summaries remain
distinct from these new smoke results. Full reference regeneration and a fresh
C++ accuracy comparison were not run for this manuscript revision.

## Figure 4 spectra (historically Figure 3)

`spectra/Work1` and `spectra/Work2` contain all eight original two-column spectral
CSVs: equivalent degree and Legendre-family power density in nT^2 per degree.
`spectrum_reanalysis.json` preserves their original file paths and hashes,
band ratios, and semilog fits over explicitly specified intervals. There are
379 rows per spectrum, on equivalent degrees 12:24:9084.

From this directory, with Python and numpy installed:

```console
python audit_spectra.py --output recomputed_spectra
```

This command was executed successfully against the packaged eight series.
It reanalyses saved spectra, without synthesizing fields or fitting a new model.
For a spectral tail, least squares fits `ln(P) = intercept + slope * degree`;
the e-folding degree is `-1/slope`. The JSON reports the actual interval,
sample size, R-squared and the number of points below the plot floor 1e-4.
The long 1812-9084 fit includes 126/115 points below that floor and gives about
202/239 degrees. Fits restricted to 1812-3900 instead give about 1153/704 degrees.
Thus no single 200-240-degree scale is assigned to the whole displayed tail in
the manuscript. A smoothly declining tail does not by itself establish a noise
level or an independently resolved wavelength.

## Regional inversion and residual audit

`application_parameters/` preserves the two small geometry and inversion setting
files. `inversion_audit.md` records the coefficient penalty, preconditioner,
stopping equations and inspected source hashes. These are format and convention
audits of the present source and saved files; no application execution manifest
or convergence log was available.

`application_residuals/` contains the saved-output audit script, aggregate JSON,
an independent stored-residual check and detailed component conventions. Running
that audit requires the original local observation and reconstruction files:

```console
python application_residuals/audit_saved_residuals.py --root PATH_TO_APPLICATIONS --output residual_check
```

Its existing statistics use all paired endpoints and no distance division.
The manuscript table reports the actual cap-frame fitted components:
Work1 (1.47, 1.74, 0.70) nT and Work2 (2.23, 1.63, 1.16) nT. Its JSON also
records geographic-frame differences, which are different horizontal quantities.
Scalar paired RMS values 0.72/1.28 nT are unchanged. These are in-sample fit
statistics, not an accuracy comparison with another model. No new original-grid
or common-physical-quantity comparison was performed.

## Building the paper

From the parent manuscript folder, run `pdflatex`, `bibtex`, `pdflatex`,
`pdflatex` on the manuscript basename. TeX Live 2026 was used for this revision.
All figures are already present under `figures/`. The Figure 5 files preserve
gold arrows and green mapped faults.

`manifest.json` records all files in this supplement except itself. The public
C++ release URL and frozen-validation/atlas identities remain in the paper's
Data and code availability statement. Full regional inversion replication still
requires the processed inputs and their preparation history; this supplement
does not provide those raw inputs or supersede their original terms.

## Geometry figure and editorial audit, 24 October 2026

The new Figure 1 and its reproducible plotting code are stored in the research
project's Geometry_Algorithm_Figure directory and included with the manuscript
source package. The four previous figures retain their filenames and appear as
Figures 2-5. editorial_audit_20261024/ records the structural/statistical audit,
figure regeneration and the new six-request Python reference smoke test with
24 controls. Those smoke checks do not replace the frozen 40,995-request C++
accuracy comparison or imply a new application inversion. The AI declaration
remains deferred at the author's explicit request and is tracked as pending.

