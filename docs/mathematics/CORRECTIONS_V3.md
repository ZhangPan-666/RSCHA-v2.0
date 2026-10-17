# Corrections in the v3 mathematical-handbook edition

## 1. Legendre cap norm

Corrected the parameter-derivative identity in Section 9.  The Dirichlet-root
formula now contains the required minus sign:

```tex
||P||_cap^2 = -(1-x0^2)/(2 nu + 1) * (partial_nu P)(partial_x P).
```

The general formula was also reordered consistently, and a `P_1` sign check
was added.

## 2. Wigner-D coefficient rotation

Under the defining identity used in Section 20, the coefficient vector obeys

```tex
c_tilde = D^T c,
```

not `D c`.  Both indexed and matrix forms now agree.

## 3. Reverse SHA--R-SCHA mapping

The former non-unique-continuation argument was replaced.  Exact harmonic
continuation to a prescribed connected domain is unique if it exists.  The
actual difficulties are non-existence for a general regional field,
ill-posedness, and non-uniqueness caused by finite or incomplete coefficient
data unless a global model class is fixed.

## 4. Boundary regularity

Section 21 now distinguishes:

- `L^2` assumptions sufficient for boundary spectral projection; and
- trace/dual-trace assumptions required for an `H^1(Omega)` weak solution.

The reconstruction and convergence statements were updated accordingly.

## 5. Handbook positioning

Prose that described the work as an article or manuscript was revised to call
it a mathematical handbook.  The LaTeX `article` document class is retained
only as a typesetting choice.
