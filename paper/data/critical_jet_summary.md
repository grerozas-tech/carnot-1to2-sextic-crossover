# Independent critical-jet re-extraction

## Method

The critical local return was re-extracted independently from the exact reduced flow. Instead of fitting sampled map values, the calculation transports a bivariate polynomial jet in the section variables `(q,p)` through the multiprecision Taylor integrator. The local event time is then solved formally, degree by degree, from `x(T_*+tau(q,p),q,p)=0`. This gives the complete return jet through total degree five in one integration.

The canonical oddification uses only the nonresonant generators required to remove the even map terms: a cubic Hamiltonian generator `G3` eliminates the quadratic map term, and a quintic generator `G5` eliminates the quartic map term. No resonant quartic Hamiltonian generator is added. This fixes the sextic gauge used in the manuscript.

## Independent linear scale

A separate multiprecision extrapolation of the detuned linear monodromy gives

- `alpha_lin = 0.0946271414514154`

The finer critical-jet transport gives

- `alpha_bal = 0.0946271414515316066063936799275`

Relative discrepancy: `1.228e-12`.

## Re-extracted normal form

Using the fixed reversible scale and the finer `hmax=0.20` transport,

- `A = 203.672358023901381`
- `B = -408.308712045131642`
- `rho = 2.001182918252685`
- `sigma = -4.733072306334130e-03`

Sixth order:

- `d = 27009.26240384307`
- `e = 271555.83498617553`
- `f = 512441.22878840013`
- `g = 75140.77917338765`
- `G = 886147.10535180639`
- `Delta6 = -289016.91057176917`

## Discretization audit

Two independent transports were recomputed with different working precisions and Taylor step bounds: `hmax=0.22` at 55 digits and `hmax=0.20` at 60 digits. Their absolute differences are

- `A`: `2.943e-09`
- `B`: `7.076e-10`
- `d`: `2.969e-06`
- `e`: `1.779e-05`
- `f`: `1.634e-05`
- `g`: `1.598e-07`
- `G`: `1.357e-06`.

The detailed comparison is stored in `critical_jet_discretization_audit.csv`.

The transformed amplitude `I_NF` is substantially more stable than the individual sextic coefficients: evaluating the two independently reconstructed canonical maps on every computed crossover point changes `I_NF` by at most `3.963e-13` relatively.

## Amplitude crossover

The reconstructed canonical chart gives the exact-return normal-form amplitude through `I_NF=(Q^2+P^2)/2`. Across `u=0` to `100`, the maximum relative deviation from the scalar amplitude curve is about `0.62%`. At the residual cycle the difference from the scalar plateau is about `0.037%`.

The first-return companions become opposite in the oddified chart to high accuracy. The relative difference between the two companion values of `I_NF` stays below `7e-6` over the full crossover range and is below `4e-8` in the inner crossover.
