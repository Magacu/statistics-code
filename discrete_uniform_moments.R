#' Moments of the Discrete Uniform Distribution U(a,b)
#' @description Compute the mean and variance for the discrete uniform distribution defined on integers a to b (inclusive).
#' @param a Integer lower bound.
#' @param b Integer upper bound (b >= a).
#' @return A list with elements `mean` and `variance`.
#' @examples
#' discrete_uniform_moments(1, 6)
#' @export
discrete_uniform_moments <- function(a, b) {
  if (any(!is.finite(c(a,b))) || (a %% 1 != 0) || (b %% 1 != 0) || b < a) {
    2
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  mean_val <- (a + b) / 2
  var_val <- ((b - a + 1)^2 - 1) / 12
  list(mean = mean_val, variance = var_val)
}