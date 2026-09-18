# Paper source and numerical reproducibility

This directory contains the complete source and numerical record for the v1.0.0 preprint.

Build:

```bash
latexmk -pdf main.tex
```

Completed audits include:

1. 120-digit semisimplicity audit of the resonant first return;
2. exact-return axial calibration of the leading quartic coefficients `A` and `rho`;
3. multiprecision residual elliptic period-two cycle, action, and Floquet angle;
4. exact-return crossover continuation from `u=0` to `u=100`;
5. elliptic action reconstruction of the diagonal sextic coefficient `G`;
6. independent multiprecision fifth-order critical-return jet transport using two discretizations;
7. canonical nonresonant oddification and nonlinear normal-form chart audit;
8. parameter-free amplitude-action and amplitude-Floquet crossover tests.

Important final correction: subleading axial action/Floquet coefficients are mixed with detuning-dependent normal-form terms at the same weighted order and are therefore **not** used to infer the frozen critical sextic tensor.

See `data/campaign_summary.md` and `data/critical_jet_final_audit.md` for the final audit summaries.
