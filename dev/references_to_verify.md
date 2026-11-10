# Methodological references and remaining provenance work

Updated 10 November 2026.

The package's GLBFP construction and Gaussian-reference L1 bandwidth heuristic
follow Carbon, M. (2024), Le fameux GLBFP, an internal working note, equations
(2)-(7) and (22)-(36). This is not presented as a peer-reviewed article or a
proof of general GLBFP asymptotic normality.

Carbon and Duchesne (2024), Multivariate frequency polygon for stationary
random fields, AISM 76, 263-287, doi:10.1007/s10463-023-00883-5, concerns the
LBFP. Carbon and Duchesne (2026; online 2025), AISM 78, 297-326,
doi:10.1007/s10463-025-00952-x, concerns LBFP asymptotic normality.
These references do not by themselves justify the general GLBFP heuristic or
confidence intervals for an arbitrary dependent-field workflow.

No separate published GLBFP methodological article is asserted without a
verified bibliographic record. The historical name compute_bi_optim() is kept
for compatibility; the current documentation states the heuristic's scope.

The bundled ashua data are retained for compatibility. Its raw hydrometric
station, exact source query, original checksum and cleaning script remain to
be established. Its provenance is not represented as fully reconstructed.
The primary real-data examples use the documented CANGRD and CDEM subsets.
