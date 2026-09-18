# Numerical campaign summary

## 1. Semisimple 1:2 resonance

Arbitrary-precision Taylor integration of the reduced flow and variational equations gives

- `mu_star = 0.6628690069305101280936908638603636487...`
- `y_star = -2.5725893074769489665217743248705277...`
- `T_quarter = 3.0179554533608596100474613978029132...`

At 120-digit working precision (Taylor order 38, `hmax=0.05`):

- `||DP + I||_F ≈ 6.24e-47`
- `|det(DP)-1| ≈ 1.89e-49`

This directly supports the semisimple identification; it is not inferred from the trace alone.

## 2. Orbital quartic calibration

Exact axial period-two actions and physical second-return Floquet exponents give

- `A_orb(delta) -> A`
- `rho_orb(delta) -> rho`

Two dyadic ladders were computed:

- ladder A: `0.020, 0.010, 0.005`
- ladder B: `0.016, 0.008, 0.004`

First Richardson values:

- `A^[1] = 203.66441, 203.66964, 203.66683, 203.67029`
- `rho^[1] = 2.0014734, 2.0012549, 2.0013681, 2.0012289`

Targets from the frozen critical jet:

- `A_jet = 203.6723580...`
- `rho_jet = 2.001182918...`

The quartic near-degeneracy is therefore visible directly in periodic-orbit observables.

## 3. Correction to the planned axial sextic reconstruction

The campaign revealed that the previous formulas attempting to reconstruct the frozen critical sextic tensor from the subleading axial action/Floquet coefficients were incomplete.

The detuned generator contains terms of the form

`K_delta = -delta I + H4 + H6 + delta H_{4,1} + gamma delta^2 I + ...`

and on axial branches `I = O(delta)`, so

- `H6 = O(delta^3)`
- `delta H_{4,1} = O(delta^3)`
- `delta^2 I = O(delta^3)`

The same terms contaminate the `delta^2` Floquet correction. Therefore the subleading axial coefficients do not isolate `d,e,f,g` unless the parameter-dependent normal form is computed to the same weighted order.

The manuscript has been corrected: axial data are now used only for the robust leading calibration of `A` and `rho`.

## 4. Residual elliptic period-two cycle

At the resonance, a representative exact `F=P^2` fixed point is

- `q = 2.568391853463491e-4`
- `p = 2.833956973721869e-3`

Its first-return companion is

- `Pq = -2.672591917458780e-4`
- `Pp = -2.703653547877403e-3`

Multiprecision observables:

- `J_exact = -1.691023541312271e-13`
- `theta_F = 1.6268450653632015e-4`
- `tr(DF) = 1.99999997353375139...`
- `|det(DF)-1| < 1e-23`

Scalar endpoint predictions:

- `J0_scalar = -1.69009171482e-13`
- `theta0_scalar = 1.62689147576e-4`

## 5. Full invariant crossover

The near-diagonal elliptic family was computed at

`u = 0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100`.

The exact normalized action and Floquet angle are compared with

- `Y_scalar = 4 R^3 - 3 R^2`
- `Z_scalar = R sqrt(2R-1)`

where `R=(1+sqrt(1+3u))/2`.

Over the full range:

- maximum normalized action deviation from the scalar curve: `< 0.81%`
- maximum normalized Floquet-angle deviation: `< 0.52%`

This directly connects the residual plateau to the outer sextic regime without fitting a crossover parameter.

## 6. Elliptic reconstruction of G

Using the exact action and the scalar equilibrium identities gives `G_ell`.
Representative values:

- `u=0.01: G_ell = 8.86023e5`
- `u=0.1:  G_ell = 8.85904e5`
- `u=1:    G_ell = 8.85186e5`
- `u=10:   G_ell = 8.82339e5`
- `u=100:  G_ell = 8.73139e5`

Frozen critical-jet target:

- `G_jet = 8.86147e5`

The inner crossover agrees at the `10^-3` level or better; the outer drift is consistent with higher-order corrections.

## 7. Independent critical-jet audit completed

The complete critical return jet through degree five has now been re-extracted independently by multiprecision bivariate jet transport with an implicit polynomial return-time correction. The canonical oddification is fixed by using only the nonresonant generators `G3` and `G5`; no resonant quartic generator is added.

Audited values:

- `alpha_lin = 0.0946271414514`
- `alpha_bal = 0.0946271414515`
- relative scale discrepancy `~1.2e-12`
- `A = 203.6723580`
- `B = -408.3087120`
- `rho = 2.001182918`
- `d = 27009.26240`
- `e = 271555.83499`
- `f = 512441.22879`
- `g = 75140.77917`
- `G = 886147.10535`

A second audit was recomputed independently with `hmax=0.22` at 55-digit working precision and compared with the finer `hmax=0.20`, 60-digit transport. The absolute differences are approximately `3.0e-6` in `d`, `1.8e-5` in `e`, `1.6e-5` in `f`, `1.6e-7` in `g`, and `1.4e-6` in `G`. The two reconstructed canonical maps change `I_NF` by at most `4.0e-13` relatively along the full crossover.

## 8. Amplitude crossover

The reconstructed nonlinear canonical chart gives `I_NF=(Q^2+P^2)/2` directly on the exact-return cycles. Across `0 <= u <= 100`, the maximum relative amplitude deviation from the scalar prediction is about `0.62%`; at the residual cycle it is about `0.037%`. Eliminating the continuation parameter gives the parameter-free relations `Y=4X^3-3X^2` and `Z=X sqrt(2X-1)`, with maximum observed departures of about `1.1%` and `0.7%`, respectively.

The only remaining extension is the full parameter-dependent normal form needed to interpret the contaminated subleading axial channel. It is not required for the present conclusions.
