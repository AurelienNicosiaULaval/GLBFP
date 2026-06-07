#' GLBFP: General Linear Blend Frequency Polygon Density Estimation
#'
#' `GLBFP` provides one-point and grid-based density estimators based on
#' ASH, LBFP and GLBFP methodology, with sparse-prefix computation,
#' leave-one-out self-support scores, visualization helpers and bandwidth
#' selection utilities.
#'
#' Main entry points:
#' - [ASH()], [LBFP()], [GLBFP()]
#' - [ASH_estimate()], [LBFP_estimate()], [GLBFP_estimate()]
#' - [compute_Di()]
#' - [compute_bi_optim()]
#'
#' Lowercase aliases such as [glbfp()] and [glbfp_estimate()] are also provided
#' for users who prefer lower-snake-case function names.
#'
#' @name GLBFP-package
"_PACKAGE"
