# Mathematical Foundations of Revised Spherical Cap Harmonic Analysis

## Current revision: November 6, 2026

This handbook derives R-SCHA on a spherical-cap shell under the boundary,
branch, normalisation, trace and convergence assumptions stated in the text.
It treats the closed-shell four-family expansion separately from classical
one-surface R-SCHA2D. Volume energy, fixed-surface gradient energy and the
Legendre-only surface approximation are distinguished.

The current build has 160 pages and 24 section modules. There are 1177 unique
labels, no unresolved references or bibliography keys, and no overfull or
underfull box warnings. These build checks do not certify mathematical
correctness. The review report records confirmed defects and their fixes,
analytical counterexamples, checks performed, and remaining scope limits.

See `notes/handbook_audit_20261106.md` for the current review. Earlier audit
reports describe earlier editions and are retained as historical records.
The old 450 numerical checks have been rerun; additional checks comprise
144 boundary/special-function comparisons, 250 energy/tensor/flux comparisons
and 588 finite-degree Wigner defining-identity comparisons. These are finite
numerical cross-checks, not a proof that the entire handbook has no errors.
New scripts and results are in `reproducibility/audit_20261106/`.

## Figures

Main basis-function and first/second derivative mode panels all use 3D cap
displays. Curve plots show all four borders. The dedicated Mehler index
diagram and repeated index-selection discussion remain removed. Numerical
illustration definitions and scripts are in `figures/`.

## Build

On the configured Windows/TeX Live installation run `./build.ps1`, or run
pdflatex, bibtex8, and two further pdflatex passes on `main.tex`. The script
uses `build/` and copies the final PDF to both root PDF filenames.

The 24 section files are assembled by `main.tex`; source filenames are not
the final printed section numbers (the one-surface section is ordered after
the three-dimensional radial theory). The assembled PDF and all coefficients
must be interpreted with the conventions and hypotheses stated in the text.

## Download and source

The current PDF is [available here](Mathematical%20Foundations%20of%20Revised%20Spherical%20Cap%20Harmonic.pdf). The complete portable LaTeX source, figures, notes and numerical check records are in [handbook-source-20261106.zip](handbook-source-20261106.zip). Extract the archive and run `build.ps1`. The v1.0 software release ZIP is a historical snapshot and contains the earlier handbook. This update does not change the numerical software or certify new software test results.
