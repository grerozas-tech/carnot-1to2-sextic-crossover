# Near-Degenerate Reversible 1:2 Resonance and Sextic Crossover

**in a (2, 3, 5, 7) Carnot Geodesic Flow**

**Author:** Gregorio Rozas Fernández  
**Affiliation:** Independent Researcher  
**Version:** v1.0.0  
**Status:** preprint release candidate, September 18, 2026  
**Repository:** https://github.com/grerozas-tech/carnot-1to2-sextic-crossover

## Overview

This repository accompanies the preprint on a near-degenerate reversible `1:2` resonance in the reduced geodesic flow of a rank-two Carnot structure with growth vector `(2,3,5,7)`.

The central mechanism is

```text
reversible 1:2 resonance
        -> quartic near-degeneracy
        -> sextic promotion in the near-diagonal sector
        -> two periodic-orbit action scales
        -> quartic-sextic crossover
        -> residual elliptic period-two cycles at resonance
```

The exact-return computations, multiprecision monodromy audit, independent critical-jet extraction, action/Floquet continuation, and machine-readable numerical tables are included for reproducibility.

## Main results

- The complete first-return monodromy at the resonant parameter is numerically consistent with the semisimple matrix `-I` at multiprecision.
- In the fixed reversible normal-form chart, the quartic Hamiltonian is extremely close to `A(Q^2-P^2)^2`.
- The resulting near-degeneracy promotes the sixth-order Hamiltonian in a distinguished near-diagonal sector.
- Axial hyperbolic period-two branches exhibit the ordinary action law `|J_H| ~ delta^2`.
- The near-diagonal elliptic branch exhibits the intermediate law `|J_E| ~ delta^(3/2)` in the sextic-dominated regime.
- Because the quartic degeneracy is only approximate, the elliptic branch crosses over to a finite-amplitude plateau.
- The exact second return contains two residual elliptic period-two cycles at resonance.
- The exact-return crossover from `u=0` to `u=100` follows the coefficient-independent scalar amplitude/action/Floquet relations at the percent or sub-percent level.
- An independent multiprecision fifth-order return-jet extraction reproduces the frozen quartic and sextic normal-form coefficients and nonlinear canonical chart.

## Important methodological correction

Subleading axial action and Floquet coefficients do **not** isolate the frozen critical sextic jet by themselves. Detuning-dependent terms enter at the same weighted order. The final analysis therefore uses axial branches only for the robust leading calibration of `A` and `rho`; the diagonal sextic combination `G` is tested through the near-diagonal elliptic crossover.

## Scientific claim boundary

This repository does **not** claim:

- a global integrability or nonintegrability theorem for the Carnot geodesic flow;
- a universal law for arbitrary reversible `1:2` resonances;
- analysis of the exactly quartic-degenerate higher-codimension limit `sigma = 0`;
- a completed theory of the associated Stokes/splitting problem.

The normalized crossover relations are coefficient-independent within the leading scalar reduction used in the paper, not generic identities of all reversible resonances.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── LICENSE-CODE
├── LICENSE-CONTENT.md
├── .zenodo.json
├── paper/
│   ├── carnot_1to2_sextic_crossover.pdf
│   ├── main.tex
│   ├── sections/
│   ├── appendices/
│   ├── figures/
│   ├── data/
│   ├── references.bib
│   ├── taylor_mp.py
│   ├── jet_transport_mp.py
│   ├── extract_critical_jet.py
│   └── ...
├── outreach/
│   ├── cuando_una_resonancia_cambia_de_escala_ES.md
│   ├── cuando_una_resonancia_cambia_de_escala_ES.pdf
│   ├── when_a_resonance_changes_scale_EN.md
│   └── when_a_resonance_changes_scale_EN.pdf
└── MANIFEST_SHA256.txt
```

## Build the paper

```bash
cd paper
latexmk -pdf main.tex
```

## Reproducibility entry points

The principal numerical scripts and audit files are in `paper/` and `paper/data/`.

- `taylor_mp.py`: arbitrary-precision Taylor integrator used in the resonance and residual-cycle audits.
- `jet_transport_mp.py`: multiprecision bivariate Taylor-jet transport.
- `extract_critical_jet.py`: critical fifth-order return-jet extraction and normal-form reconstruction.
- `numerics_double.py`: exact local return, event-corrected monodromy, continuation, and section-action quadrature.
- `paper/data/campaign_summary.md`: numerical campaign summary.
- `paper/data/critical_jet_final_audit.md`: final independent critical-jet and nonlinear-chart audit.

## Popular-language explanation

Nontechnical explanations are included in Spanish and English under `outreach/`.

In one sentence:

> A near-cancellation of the leading nonlinear restoring term makes the next order dynamically dominant, producing two periodic-orbit scales and a measurable crossover.

## Licensing

- Source code and numerical scripts: MIT.
- Original paper text, explanatory text, original figures, and original numerical presentation: CC BY 4.0.

## Citation and DOI

A `CITATION.cff` file is included. The version DOI will be added here after the v1.0.0 Zenodo deposit is published.
