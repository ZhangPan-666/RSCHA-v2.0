# Validation report

## Deliverable

- Title: *Mathematical Foundations of Revised Spherical Cap Harmonic Analysis*
- Edition: corrected mathematical-handbook edition
- Date: July 31, 2026
- PDF pages: 166
- Main section modules: 24
- Section-source lines: 17,846
- LaTeX labels: 1,207, all unique
- Cross-reference occurrences: 237, addressing 292 label targets
- Bibliography database entries: 12
- Citation occurrences: 14, using 8 bibliography keys

## Corrections made in this edition

1. Corrected the parameter-derivative formula for the associated-Legendre cap
   norm.  Both the general diagonal-limit identity and the Dirichlet-root
   simplification now have the correct sign and term order.
2. Added a direct sign check using `m=0`, `x0=0`, and `nu=1`, for which the
   corrected formula gives `1/3`, equal to the direct integral of `P_1^2`.
3. Reconciled the Wigner-`D` defining identity with the coefficient-space
   transformation.  Under the declared index convention, column coefficients
   transform as `c_tilde = D^T c`.
4. Replaced the incorrect claim of non-unique exact harmonic continuation.
   The handbook now states that continuation to a prescribed connected domain
   is unique if it exists, while the actual reverse-mapping obstacles are
   non-existence, severe ill-posedness, and incomplete finite-dimensional data.
5. Distinguished boundary `L^2` spectral projection from the trace and
   dual-trace hypotheses required for an `H^1(Omega)` weak solution.
6. Updated the reconstruction and convergence statements so that `L^2`
   boundary reconstruction is not automatically identified with finite-energy
   convergence.

## Previously audited corrections retained

- Removal of Pan Zhang citations and attribution.
- Concise placement of classical R-SCHA2D after the three-dimensional
  eigenfunction construction.
- Corrected `V2` zero-mode projection normalisation.
- Corrected Mehler bilinear Wronskian sign and order.
- Corrected continuous fixed-surface full-gradient Gram-form null space.
- Quotient trace-space and Lions--Magenes-type weak Neumann formulation.
- Removal of the former coordinate-frame illustration while retaining the
  analytical coordinate formulas and cap-rotation figure.

## Structural checks

- Duplicate labels: 0
- Missing cross-reference labels after final compilation: 0
- Missing bibliography keys: 0
- Missing referenced image files: 0
- LaTeX undefined-reference or citation warnings: 0
- LaTeX overfull/underfull box warnings: 0
- PDF openable: yes
- PDF encrypted: no
- Searchable text: yes
- Page size: US Letter

## Mathematical checks

- The corrected Legendre cap-norm formula was checked at
  `m=0`, `x0=0`, `nu=1`.  Numerical parameter differentiation gives
  `partial_nu P_nu(0)|_{nu=1}=-1` and `partial_x P_1(0)=1`, so the corrected
  negative-sign formula gives `1/3`, agreeing with direct quadrature.
- The Wigner-`D` coefficient transformation was re-derived by substituting the
  defining basis identity into a scalar degree block and collecting the local
  basis index.  This gives `c_tilde_q = sum_m D_mq c_m`, equivalently
  `c_tilde = D^T c` for column storage.
- The reverse-mapping discussion was checked against real analyticity and
  unique continuation of harmonic functions on connected domains.
- The boundary-data discussion was checked for consistency with the weak
  formulation and data spaces already declared in Section 4.

## PDF rendering checks

- Poppler rendering completed successfully for all 166 pages at 110 dpi.
- The title page and every modified passage were additionally rendered at
  180 dpi and visually inspected.
- Inspected modification pages include the Legendre cap norm, Wigner-`D`
  rotation, reverse mapping, boundary-data spaces, reconstruction conditions,
  and convergence hierarchy.
- No clipping, overlap, broken glyphs, or malformed equation numbering was
  observed on the inspected pages.

## Mathematical limits retained explicitly

This edition does not convert open problems into theorems.  The completeness
or frame status of the commonly used triangular Mehler index family, sharp
mixed-edge regularity, and a full functional-analytic basis theorem for the
combined classical R-SCHA2D family remain qualified in the handbook.
