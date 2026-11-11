# Stationary-mixing LBFP intervals: numerical checks

The implementation follows Carbon and Duchesne (2026), *Asymptotic normality
of multivariate frequency polygons for stationary random fields*, Annals of
the Institute of Statistical Mathematics 78, 297-326,
<https://doi.org/10.1007/s10463-025-00952-x>. The formula and hypotheses were
checked against the authors' revised manuscript R2 dated 3 June 2025, obtained
from the authors' correspondence. Equation numbers below refer to that
manuscript; the version-of-record text has not been compared in full.

For the LBFP, manuscript equation (13) estimates the scaled asymptotic variance
by `fhat(x) * sum_j(c_j(x)^2)`. The weights are multilinear products, hence the
sum of their squares equals `product_k((1-u_k)^2 + u_k^2)`. Equation (16)
divides by `n * product(b)` before taking a square root. The package uses the
interpolation phase on its counting grid, with unoccupied cells completed by
zero. The independent tests construct histogram counts and corner weights
directly, using all three numerical engines, dimensions N=1,2,3 and d=1,2.

The theorem concerns a stationary, mixing random field on expanding rectangular
regions with bounded aspect ratios, regular joint and marginal densities,
positive density at the fixed evaluation point, and specified rates for mixing
and bandwidths. The corollary and theorem leading to the displayed interval
add a negligible-bias requirement. The code checks the geometry and explicit
parameters; it cannot establish the probabilistic or asymptotic assumptions.
A zero plug-in density gives `NA`, not a zero-width interval.

The manuscript also provides a heuristic variance expression for GLBFP with
m>1 and discusses simulations. Its general asymptotic justification is left
open by its authors. It is not implemented as a theoretically validated
interval here. ASH is also outside this implementation.

## Simulation protocol

`simulate-mixing-fields.R` uses 500 replications per scenario and N, with seed
20261111. Observations have marginal standard normal distribution. The query
is x=0.15, bandwidth b=0.3, anchor -5, and the nominal level is 95%. The
rectangles are length 900, 30 by 30, and 10 by 10 by 10 for N=1,2,3. Scenarios
are independent observations, a finite-range moving average with coefficient
0.5, and a separable stationary Gaussian autoregressive construction with
correlation parameter 0.9 along each coordinate. The latter has strong
finite-sample correlation. For each simulated field the script records the
estimate, interval, true marginal density and coverage; non-finite outputs
stop the experiment rather than being silently discarded.

| N | Scenario | n | Coverage | Monte Carlo standard error |
| ---: | --- | ---: | ---: | ---: |
| 1 | independent | 900 | 96.8% | 0.79 percentage point |
| 2 | independent | 900 | 97.2% | 0.74 percentage point |
| 3 | independent | 1000 | 96.2% | 0.86 percentage point |
| 1 | finite-range MA | 900 | 94.4% | 1.03 percentage points |
| 2 | finite-range MA | 900 | 97.0% | 0.76 percentage point |
| 3 | finite-range MA | 1000 | 95.6% | 0.92 percentage point |
| 1 | strong AR | 900 | 85.0% | 1.60 percentage points |
| 2 | strong AR | 900 | 47.2% | 2.23 percentage points |
| 3 | strong AR | 1000 | 29.2% | 2.03 percentage points |

These results reveal substantial finite-sample undercoverage under strong
correlation. They do not disprove an asymptotic theorem. In particular, these
runs use a single fixed bandwidth and do not verify a bandwidth sequence or
any asymptotic rate condition. Finite-sample validity is not claimed.

`check-gaussian-generator.R` separately compares simulated covariance matrices
of small separable AR fields against the theoretical values over 5000 draws.
The maximum absolute discrepancies by N are recorded in `generator-check.csv`.
The CSV files contain all elementary replications and summary statistics;
`session-info.txt` records the computing environment.
