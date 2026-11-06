# Audited revision notes

This edition was rebuilt from the detailed pure-mathematics project.

## Removed

- Every Pan Zhang citation, bibliography entry, and dissertation attribution.
- The former full Part that treated R-SCHA2D as a collapsed-shell limit.
- The low-quality geocentric/spherical/local coordinate-frame illustration.
- Residual numerical implementation and stability prescriptions.

## Rewritten

- Classical R-SCHA2D is now a concise standalone branch following the
  three-dimensional construction.  It uses the lateral derivative root
  condition of Thébault (2008), the two Euler radial branches for `V1`, and
  the separate repeated degree `-1/2` continuation for `V2`.
- The summary, introduction, abstract, and notation register now distinguish
  classical one-surface R-SCHA2D from the closed-shell three-dimensional
  R-SCHA spectrum.

## Corrected mathematical statements

- `V2` constant-mode boundary projection normalisation.
- Mehler bilinear Wronskian sign.
- Continuous fixed-surface Gram-form null space.
- Mixed-boundary trace spaces and weak Neumann pairing.
- Scope of the infinite `V2` series as an `H^1` weak solution.
- Non-vanishing of the lateral Mehler factor appearing in coefficient
  denominators.

## Retained typography

- US Letter paper with one-inch margins.
- 11 pt Computer Modern text.
- Enlarged, vertically arranged Legendre and Mehler harmonic panels.
- Article-style contents, headings, captions, and bibliography.

## Corrected after mathematical review

- Reversed the sign and term order in the parameter-derivative formula for
  the Legendre cap norm; added a direct $P_1$ sign check.
- Corrected the Wigner-$D$ coefficient rotation to
  `c_tilde = D^T c` under the handbook's defining identity.
- Replaced the claim of non-unique exact harmonic continuation by the correct
  distinction: continuation may fail to exist, is unique when it exists on a
  prescribed connected domain, and is generally severely ill posed.
- Distinguished $L^2$ boundary spectral projection from the
  $H^{1/2}$ / dual-trace hypotheses required for an $H^1$ weak solution.
