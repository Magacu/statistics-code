#' Statistical Inference Library
#'
#' A collection of functions for statistical inference, including confidence intervals
#' and hypothesis testing for means, proportions, and variances, as well as
#' non-parametric tests like Chi-square goodness of fit and independence.
#' Based on standard statistical inference curriculum.
#'
#' @name inference_library
NULL

#' Confidence Interval for the Mean
#' @family confidence-intervals
#' @description Calculates the confidence interval for a population mean μ.
#' Supports scenarios with known or unknown population variance, and large or small sample sizes.
#'
#' @details This function implements the formulas for CI of the mean based on the
#' central limit theorem and t-distributions.
#' \itemize{
#'   \item Known variance: Uses Normal distribution (z-score).
#'   \item Unknown variance, large n (>30): Uses Normal approximation (z-score).
#'   \item Unknown variance, small n (<=30): Uses Student's t-distribution.
#' }
#'
#' @param x Numeric vector of data. If NULL, summary statistics (mean_x, n, sd_x) must be provided.
#' @param conf_level Numeric value between 0 and 1 indicating the confidence level (1 - α). Default is 0.95.
#' @param sigma_pop Numeric (optional). The known population standard deviation (σ). If NULL, sample standard deviation is used.
#' @param mean_x Numeric (optional). Sample mean (if x is not provided).
#' @param n Integer (optional). Sample size (if x is not provided).
#' @param sd_x Numeric (optional). Sample standard deviation (if x is not provided and sigma_pop is NULL).
#' @return A named numeric vector containing the lower and upper bounds of the confidence interval.
#'
#' @examples
#' # From raw data
#' data <- rnorm(25, mean = 10, sd = 2)
#' ci_mean(data, conf_level = 0.95)
#'
#' # From summary statistics (n=144, mean=160, var=100 -> sd=10) - Example from PDF Page 34
#' ci_mean(mean_x = 160, n = 144, sigma_pop = 10, conf_level = 0.95)
#'
#' @export
ci_mean <- function(x = NULL, conf_level = 0.95, sigma_pop = NULL, mean_x = NULL, n = NULL, sd_x = NULL) {
  if (!is.null(x)) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
  }
  
  if (is.null(n) || is.null(mean_x)) stop("Must provide raw data 'x' or summary stats 'n' and 'mean_x'.")
  
  alpha <- 1 - conf_level
  
  # Logic for selecting distribution based on PDF rules
  if (!is.null(sigma_pop)) {
    # Known variance -> Normal Distribution
    crit_val <- qnorm(1 - alpha / 2)
    margin_error <- crit_val * (sigma_pop / sqrt(n))
    method <- "Normal (Known Variance)"
  } else {
    # Unknown variance
    if (is.null(sd_x)) stop("Sample standard deviation 'sd_x' required if 'sigma_pop' is unknown.")
    
    if (n > 30) {
      # Large sample -> Normal Approximation
      crit_val <- qnorm(1 - alpha / 2)
      margin_error <- crit_val * (sd_x / sqrt(n))
      method <- "Normal Approx (Large n)"
    } else {
      # Small sample -> T-Student
      crit_val <- qt(1 - alpha / 2, df = n - 1)
      margin_error <- crit_val * (sd_x / sqrt(n))
      method <- "T-Student (Small n)"
    }
  }
  
  ci <- c(lower = mean_x - margin_error, upper = mean_x + margin_error)
  attr(ci, "method") <- method
  return(ci)
}

#' Confidence Interval for Variance
#' @family confidence-intervals
#' @description Calculates the confidence interval for a population variance σ².
#' Uses the Chi-square distribution.
#'
#' @details Assumes the underlying population is normally distributed.
#' Formula: \eqn{[ \frac{(n-1)S^2}{\chi^2_{1-\alpha/2}}, \frac{(n-1)S^2}{\chi^2_{\alpha/2}} ]}
#'
#' @param x Numeric vector of data. If NULL, summary statistics (var_x, n) must be provided.
#' @param conf_level Numeric value between 0 and 1. Default is 0.95.
#' @param var_x Numeric (optional). Sample variance (S²) (quasi-variance).
#' @param n Integer (optional). Sample size.
#' @return A named numeric vector with lower and upper bounds for σ².
#'
#' @examples
#' # Example from PDF Page 18: variance of salaries
#' # n=25, sd=1200 -> var = 1440000
#' ci_variance(n = 25, var_x = 1440000, conf_level = 0.95)
#'
#' @export
ci_variance <- function(x = NULL, conf_level = 0.95, var_x = NULL, n = NULL) {
  if (!is.null(x)) {
    n <- length(x)
    var_x <- var(x)
  }
  
  if (is.null(n) || is.null(var_x)) stop("Must provide data 'x' or 'n' and 'var_x'.")
  
  alpha <- 1 - conf_level
  df <- n - 1
  
  chi_upper <- qchisq(1 - alpha / 2, df)
  chi_lower <- qchisq(alpha / 2, df)
  
  lower <- (df * var_x) / chi_upper
  upper <- (df * var_x) / chi_lower
  
  return(c(lower = lower, upper = upper))
}

#' Confidence Interval for Proportion
#' @family confidence-intervals
#' @description Calculates the confidence interval for a population proportion p.
#' Uses the Normal approximation for large samples.
#'
#' @details
#' Formula: \eqn{\hat{p} \pm z_{\alpha/2} \sqrt{\frac{\hat{p}(1-\hat{p})}{n}}}
#'
#' @param x Numeric vector of binary data (0/1) OR integer count of successes.
#' @param n Integer. Sample size (required if x is a count).
#' @param conf_level Numeric value between 0 and 1. Default is 0.95.
#' @return A named numeric vector with lower and upper bounds.
#'
#' @examples
#' # Example from PDF Page 25: 27 infested fruits out of 150
#' ci_proportion(x = 27, n = 150, conf_level = 0.95)
#'
#' @export
ci_proportion <- function(x, n = NULL, conf_level = 0.95) {
  if (length(x) > 1) {
    # x is raw data vector
    n <- length(x)
    p_hat <- mean(x) # Proportion of 1s
  } else {
    # x is number of successes
    if (is.null(n)) stop("If 'x' is a count, 'n' must be provided.")
    p_hat <- x / n
  }
  
  alpha <- 1 - conf_level
  z_score <- qnorm(1 - alpha / 2)
  
  margin_error <- z_score * sqrt((p_hat * (1 - p_hat)) / n)
  
  return(c(lower = max(0, p_hat - margin_error), upper = min(1, p_hat + margin_error)))
}

#' Confidence Interval for Difference of Means
#' @family confidence-intervals
#' @description Calculates CI for the difference between two population means (μ1 - μ2).
#' Handles equal/unequal variances and large/small samples.
#'
#' @details
#' \itemize{
#'   \item Large samples: Uses Normal approximation.
#'   \item Small samples, equal variances: Uses pooled variance and T-distribution.
#'   \item Small samples, unequal variances: Uses Welch's approximation for degrees of freedom.
#' }
#'
#' @param x1 Numeric vector for sample 1 or list of summary stats (mean, sd, n).
#' @param x2 Numeric vector for sample 2 or list of summary stats.
#' @param conf_level Numeric. Default 0.95.
#' @param var_equal Logical. Assume equal variances? Default FALSE.
#' @return A named numeric vector with bounds for the difference.
#'
#' @examples
#' # Example data
#' set.seed(123)
#' s1 <- rnorm(35, 10, 2)
#' s2 <- rnorm(40, 12, 2.5)
#' ci_diff_means(s1, s2, var_equal = FALSE)
#'
#' @export
ci_diff_means <- function(x1, x2, conf_level = 0.95, var_equal = FALSE) {
  # Handle input (vectors vs summary stats)
  if (is.numeric(x1) && length(x1) > 1) {
    n1 <- length(x1); m1 <- mean(x1); s1 <- sd(x1)
  } else {
    n1 <- x1$n; m1 <- x1$mean; s1 <- x1$sd
  }
  
  if (is.numeric(x2) && length(x2) > 1) {
    n2 <- length(x2); m2 <- mean(x2); s2 <- sd(x2)
  } else {
    n2 <- x2$n; m2 <- x2$mean; s2 <- x2$sd
  }
  
  alpha <- 1 - conf_level
  diff_means <- m1 - m2
  
  # Logic from PDF Page 21
  if (n1 + n2 > 60) { # Threshold for "Large" is typically 30 per group or sum > 30-60
    crit_val <- qnorm(1 - alpha / 2)
    se <- sqrt(s1^2/n1 + s2^2/n2)
  } else {
    # Small samples
    if (var_equal) {
      df <- n1 + n2 - 2
      sp_sq <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / df
      se <- sqrt(sp_sq * (1/n1 + 1/n2))
      crit_val <- qt(1 - alpha / 2, df)
    } else {
      # Welch
      num <- (s1^2/n1 + s2^2/n2)^2
      den <- ((s1^2/n1)^2 / (n1 - 1)) + ((s2^2/n2)^2 / (n2 - 1))
      df <- num / den
      se <- sqrt(s1^2/n1 + s2^2/n2)
      crit_val <- qt(1 - alpha / 2, df)
    }
  }
  
  return(c(lower = diff_means - crit_val * se, upper = diff_means + crit_val * se))
}

#' Hypothesis Test for the Mean (One Sample)
#' @family hypothesis-testing
#' @description Performs a hypothesis test for a population mean μ.
#'
#' @details Computes the test statistic and p-value for H0: μ = mu_0.
#' Supports one-sided ("less", "greater") and two-sided ("two.sided") alternatives.
#'
#' @param x Numeric vector of data.
#' @param mu0 Numeric. The hypothesized mean value.
#' @param sigma_pop Numeric (optional). Known population standard deviation.
#' @param alternative Character string. "two.sided", "less", or "greater".
#' @return A list containing statistic, p.value, and method.
#'
#' @examples
#' data <- c(2.2, 2.66, 2.74, 3.41, 2.46) # Subset of PDF page 26
#' test_mean(data, mu0 = 2, alternative = "greater")
#'
#' @export
test_mean <- function(x, mu0, sigma_pop = NULL, alternative = "two.sided") {
  n <- length(x)
  mean_x <- mean(x)
  
  if (!is.null(sigma_pop)) {
    stat <- (mean_x - mu0) / (sigma_pop / sqrt(n))
    dist <- "norm"
    method <- "Z-test (Known Variance)"
  } else {
    sd_x <- sd(x)
    stat <- (mean_x - mu0) / (sd_x / sqrt(n))
    dist <- "t"
    method <- "T-test (Unknown Variance)"
  }
  
  if (dist == "norm") {
    if (alternative == "less") p_val <- pnorm(stat)
    else if (alternative == "greater") p_val <- 1 - pnorm(stat)
    else p_val <- 2 * (1 - pnorm(abs(stat)))
  } else {
    if (alternative == "less") p_val <- pt(stat, df = n - 1)
    else if (alternative == "greater") p_val <- 1 - pt(stat, df = n - 1)
    else p_val <- 2 * (1 - pt(abs(stat), df = n - 1))
  }
  
  return(list(statistic = stat, p.value = p_val, method = method, alternative = alternative))
}

#' Chi-Square Goodness of Fit Test
#' @family hypothesis-testing
#' @description Tests whether observed frequency data fits a theoretical distribution.
#'
#' @details Corresponds to the "Bondad de Ajuste" section in the PDF.
#' Statistic: \eqn{\chi^2 = \sum \frac{(O_i - E_i)^2}{E_i}}
#'
#' @param observed Numeric vector of observed frequencies.
#' @param expected Numeric vector of expected probabilities or frequencies.
#' @return A list with statistic, degrees of freedom, and p-value.
#'
#' @examples
#' # Example PDF Page 69 (Geiger counter)
#' obs <- c(200, 220, 150, 68, 25, 10, 4)
#' # Need to calculate expected probabilities for Poisson...
#' # This function serves as the generic wrapper for chisq.test
#' test_goodness_fit(obs, p = rep(1/7, 7)) # Dummy expected
#'
#' @export
test_goodness_fit <- function(observed, expected) {
  # Normalize expected if they are probabilities
  if (sum(expected) > 1.01 || sum(expected) < 0.99) {
    # Assume frequencies, so normalize to probs
    expected <- expected / sum(expected)
  }
  
  res <- chisq.test(x = observed, p = expected)
  return(list(statistic = res$statistic, df = res$parameter, p.value = res$p.value))
}

#' Chi-Square Test of Independence
#' @family hypothesis-testing
#' @description Tests for independence between two categorical variables using a contingency table.
#'
#' @details Corresponds to the "Contingency Tables" section in the PDF.
#'
#' @param table_matrix A numeric matrix or table representing the contingency table (Observed frequencies).
#' @return A list with statistic, df, p-value, and expected matrix.
#'
#' @examples
#' # Example PDF Page 74: Pear Rust vs Scab
#' data <- matrix(c(198, 28, 62, 39, 6, 12, 105, 15, 35), nrow=3, byrow=TRUE)
#' test_independence(data)
#'
#' @export
test_independence <- function(table_matrix) {
  res <- chisq.test(table_matrix)
  return(list(
    statistic = res$statistic,
    df = res$parameter,
    p.value = res$p.value,
    expected = res$expected,
    observed = res$observed
  ))
}

#' Hypothesis Test for Equality of Variances (F-test)
#' @family hypothesis-testing
#' @description Tests if two population variances are equal.
#'
#' @details Corresponds to the test for \eqn{\sigma_1^2 = \sigma_2^2}.
#' Statistic: \eqn{F = S_1^2 / S_2^2}
#'
#' @param x1 Numeric vector sample 1.
#' @param x2 Numeric vector sample 2.
#' @param alternative Character string. "two.sided", "less", "greater".
#' @return List with F-statistic, df, and p-value.
#'
#' @export
test_var_ratio <- function(x1, x2, alternative = "two.sided") {
  res <- var.test(x1, x2, alternative = alternative)
  return(list(statistic = res$statistic, df = res$parameter, p.value = res$p.value, conf.int = res$conf.int))
}