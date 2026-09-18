# Final critical-jet audit

The critical first-return jet was reconstructed directly from the reduced flow by multiprecision bivariate Taylor-jet transport through total section degree five. The local return time is solved as a polynomial from the event equation `x(T_*+tau(q,p),q,p)=0`.

Two independent transport settings were compared:

- audit A: 55-digit working precision, Taylor time order 22, `hmax=0.22`;
- audit B: 60-digit working precision, Taylor time order 22, `hmax=0.20`.

The linear part is subsequently fixed to the independently audited semisimple value `-I`; the pre-forcing linear errors are about `7.0e-11` and `5.4e-12`, respectively.

The nonlinear gauge is fixed by canonical nonresonant oddification only: a cubic Hamiltonian generator removes the quadratic map term and a quintic Hamiltonian generator removes the quartic map term. No resonant quartic generator is added.

Using the finer transport and the independently calibrated reversible scale `alpha_lin=0.0946271414514154`, the frozen coefficients are

- `A = 203.67235802390138`
- `B = -408.30871204513164`
- `rho = 2.0011829182526855`
- `sigma = -0.004733072306335014`
- `d = 27009.26240384307`
- `e = 271555.8349861755`
- `f = 512441.2287884001`
- `g = 75140.77917338767`
- `G = 886147.1053518064`
- `Delta6 = -289016.9105717692`

The full coefficient-by-coefficient discretization comparison is stored in `critical_jet_discretization_audit.csv`. The differences between the two transports are at or below about `1e-10` relative for the sextic coefficients and substantially smaller for the quartic coefficients.

The nonlinear canonical maps reconstructed from the two transports were applied independently to every exact-return cycle on the crossover grid `0 <= u <= 100`. The resulting normal-form amplitude `I_NF=(Q^2+P^2)/2` differs by at most `3.97e-13` relatively between the two maps.

This closes the independent critical-jet, canonical oddification, and nonlinear-chart audit required by the present manuscript. The only optional extension is a parameter-dependent normal form sufficient to disentangle detuning-dependent weighted-order terms from the frozen sextic tensor on the axial branches.
