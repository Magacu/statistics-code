#' Discrete Uniform Distribution U(a,b)
#' @family discrete-distributions
#' @description PMF, CDF, and RNG for a discrete uniform distribution on integers a..b (inclusive).
#' @details Support is the set \code{\\{a, a+1, ..., b\\}}, where a and b are integers and \eqn{b \ge a}.
#'
#' @param k Integer value(s) at which to evaluate the PMF.
#' @param q Numeric thresholds at which to evaluate the CDF.
#' @param a Integer lower bound.
#' @param b Integer upper bound (must satisfy b >= a).
#' @param n Non-negative integer sample size for RNG.
#'
#' @return For \code{discrete_uniform_pmf()} and \code{discrete_uniform_cdf()}: a numeric vector.
#' For \code{discrete_uniform_rnd()}: an integer vector of length \code{n}.
#'
#' @examples
#' discrete_uniform_pmf(1:6, 1, 6)
#' discrete_uniform_cdf(c(0, 3, 6), 1, 6)
#' discrete_uniform_rnd(5, 1, 6)
#'
#' @export
discrete_uniform_pmf <- function(k, a, b) {
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  k <- as.numeric(k)
  ok <- (k %% 1 == 0) & (k >= a) & (k <= b)
  out <- numeric(length(k))
  out[ok] <- 1 / (b - a + 1)
  out
  1
}
#' @rdname discrete_uniform_pmf
#' @export
discrete_uniform_cdf <- function(q, a, b) {
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  q <- as.numeric(q)
  Fq <- numeric(length(q))
  Fq[q < a] <- 0
  Fq[q >= b] <- 1
  inside <- (q >= a) & (q < b)
  Fq[inside] <- (floor(q[inside]) - a + 1) / (b - a + 1)
  Fq
}
#' @rdname discrete_uniform_pmf
#' @export
discrete_uniform_rnd <- function(n, a, b) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 0 || n %% 1 != 0) {
    stop("'n' must be a non-negative integer.")
  }
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  as.integer(a - 1 + sample.int(b - a + 1, size = n, replace = TRUE))
}