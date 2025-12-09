#' @param x Vector de datos. Si es NULL, usar estadísticos (mean_x, n, sd_x).
#' @param conf_level Nivel de confianza (0-1). Default 0.95.
#' @param sigma_pop Desviación estándar poblacional conocida (opcional).
#' @param mean_x Media muestral (si x es NULL).
#' @param n Tamaño de muestra (si x es NULL).
#' @param sd_x Desviación estándar muestral (si x es NULL).
#' @export
ci_mean <- function(x = NULL, conf_level = 0.95, sigma_pop = NULL, mean_x = NULL, n = NULL, sd_x = NULL) {
  if (!is.null(x)) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
  }
  
  if (is.null(n) || is.null(mean_x)) stop("Provide raw 'x' or summary stats 'n' and 'mean_x'.")
  
  alpha <- 1 - conf_level
  
  if (!is.null(sigma_pop)) {
    # Varianza conocida (Normal)
    crit_val <- qnorm(1 - alpha / 2)
    margin_error <- crit_val * (sigma_pop / sqrt(n))
    method <- "Normal (Known Variance)"
  } else {
    # Varianza desconocida
    if (is.null(sd_x)) stop("Sample SD 'sd_x' required if 'sigma_pop' is unknown.")
    
    if (n > 30) {
      crit_val <- qnorm(1 - alpha / 2)
      margin_error <- crit_val * (sd_x / sqrt(n))
      method <- "Normal Approx (Large n)"
    } else {
      crit_val <- qt(1 - alpha / 2, df = n - 1)
      margin_error <- crit_val * (sd_x / sqrt(n))
      method <- "T-Student (Small n)"
    }
  }
  
  ci <- c(lower = mean_x - margin_error, upper = mean_x + margin_error)
  attr(ci, "method") <- method
  return(ci)
}

#' @param x Vector de datos. Si es NULL, usar (var_x, n).
#' @param conf_level Nivel de confianza. Default 0.95.
#' @param var_x Varianza muestral (S^2).
#' @param n Tamaño de muestra.
#' @export
ci_variance <- function(x = NULL, conf_level = 0.95, var_x = NULL, n = NULL) {
  if (!is.null(x)) {
    n <- length(x)
    var_x <- var(x)
  }
  
  if (is.null(n) || is.null(var_x)) stop("Provide data 'x' or 'n' and 'var_x'.")
  
  alpha <- 1 - conf_level
  df <- n - 1
  
  chi_upper <- qchisq(1 - alpha / 2, df)
  chi_lower <- qchisq(alpha / 2, df)
  
  lower <- (df * var_x) / chi_upper
  upper <- (df * var_x) / chi_lower
  
  return(c(lower = lower, upper = upper))
}

#' @param x Vector binario o conteo de éxitos.
#' @param n Tamaño de muestra (requerido si x es conteo).
#' @param conf_level Nivel de confianza. Default 0.95.
#' @export
ci_proportion <- function(x, n = NULL, conf_level = 0.95) {
  if (length(x) > 1) {
    n <- length(x)
    x_count <- sum(x)
    p_hat <- mean(x)
  } else {
    if (is.null(n)) stop("If 'x' is a count, 'n' must be provided.")
    x_count <- x
    p_hat <- x_count / n
  }
  
  method <- ""
  
  if (n > 30) {
    alpha <- 1 - conf_level
    z_score <- qnorm(1 - alpha / 2)
    margin_error <- z_score * sqrt((p_hat * (1 - p_hat)) / n)
    lower <- max(0, p_hat - margin_error)
    upper <- min(1, p_hat + margin_error)
    ci <- c(lower = lower, upper = upper)
    method <- "Normal Approx (Large n)"
  } else {
    res <- binom.test(x_count, n, conf.level = conf_level)
    ci <- res$conf.int
    names(ci) <- c("lower", "upper")
    method <- "Exact Binomial (Small n)"
  }
  
  attr(ci, "method") <- method
  return(ci)
}

#' @param x1 Muestra 1 (vector) o lista de estadísticos (mean, sd, n).
#' @param x2 Muestra 2 (vector) o lista de estadísticos.
#' @param conf_level Nivel de confianza. Default 0.95.
#' @param var_equal Asumir varianzas iguales (TRUE/FALSE).
#' @param sigma1 SD poblacional grupo 1 (opcional).
#' @param sigma2 SD poblacional grupo 2 (opcional).
#' @export
ci_diff_means <- function(x1, x2, conf_level = 0.95, var_equal = FALSE, sigma1 = NULL, sigma2 = NULL) {
  # Parse input x1
  if (is.numeric(x1) && length(x1) > 1) {
    n1 <- length(x1); m1 <- mean(x1); s1 <- sd(x1)
  } else {
    n1 <- x1$n; m1 <- x1$mean; s1 <- x1$sd
  }
  
  # Parse input x2
  if (is.numeric(x2) && length(x2) > 1) {
    n2 <- length(x2); m2 <- mean(x2); s2 <- sd(x2)
  } else {
    n2 <- x2$n; m2 <- x2$mean; s2 <- x2$sd
  }
  
  alpha <- 1 - conf_level
  diff_means <- m1 - m2
  method <- ""
  
  if (!is.null(sigma1) && !is.null(sigma2)) {
    crit_val <- qnorm(1 - alpha / 2)
    se <- sqrt(sigma1^2/n1 + sigma2^2/n2)
    method <- "Normal (Known Variances)"
  } else {
    if (n1 + n2 > 30) { 
      crit_val <- qnorm(1 - alpha / 2)
      se <- sqrt(s1^2/n1 + s2^2/n2)
      method <- "Normal Approx (Large n, Unknown Var)"
    } else {
      if (var_equal) {
        df <- n1 + n2 - 2
        sp_sq <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / df
        se <- sqrt(sp_sq * (1/n1 + 1/n2))
        crit_val <- qt(1 - alpha / 2, df)
        method <- "T-Student (Small n, Pooled Var)"
      } else {
        num <- (s1^2/n1 + s2^2/n2)^2
        den <- ((s1^2/n1)^2 / (n1 - 1)) + ((s2^2/n2)^2 / (n2 - 1))
        df <- num / den
        se <- sqrt(s1^2/n1 + s2^2/n2)
        crit_val <- qt(1 - alpha / 2, df)
        method <- "T-Student (Small n, Welch)"
      }
    }
  }
  
  ci <- c(lower = diff_means - crit_val * se, upper = diff_means + crit_val * se)
  attr(ci, "method") <- method
  return(ci)
}

#' @param x Vector de datos.
#' @param mu0 Media hipotética.
#' @param sigma_pop Desviación estándar conocida (opcional).
#' @param alternative "two.sided", "less", o "greater".
#' @export
test_mean <- function(x = NULL, mu0, sigma_pop = NULL, alternative = "two.sided", mean_x = NULL, n = NULL, sd_x = NULL) {
  if (!is.null(x)) {
    n <- length(x)
    mean_x <- mean(x)
    sd_x <- sd(x)
  }
  
  if (is.null(n) || is.null(mean_x)) stop("Must provide raw data 'x' or summary stats 'n' and 'mean_x'.")
  
  # Logic based on PDF Page 49
  if (!is.null(sigma_pop)) {
    # Known Variance -> Z-test
    stat <- (mean_x - mu0) / (sigma_pop / sqrt(n))
    dist <- "norm"
    method <- "Z-test (Known Variance)"
  } else {
    # Unknown Variance
    if (is.null(sd_x)) stop("Must provide 'sd_x' if variance is unknown.")
    
    if (n > 30) {
      # Large Sample -> Z-test approx
      stat <- (mean_x - mu0) / (sd_x / sqrt(n))
      dist <- "norm"
      method <- "Z-test Approx (Unknown Variance, Large n)"
    } else {
      # Small Sample -> T-test
      stat <- (mean_x - mu0) / (sd_x / sqrt(n))
      dist <- "t"
      method <- "T-test (Unknown Variance, Small n)"
    }
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


#' @param x Conteo de éxitos o vector.
#' @param n Tamaño de muestra (requerido si x es conteo).
#' @param p0 Proporción hipotética.
#' @param alternative "two.sided", "less", "greater".
#' @export
test_proportion <- function(x, n = NULL, p0 = 0.5, alternative = "two.sided") {
  if (length(x) > 1) {
    n <- length(x)
    x_count <- sum(x)
    p_hat <- mean(x)
  } else {
    if (is.null(n)) stop("If 'x' is a count, 'n' must be provided.")
    x_count <- x
    p_hat <- x_count / n
  }
  
  if (n > 30) {
    se <- sqrt(p_hat * (1 - p_hat) / n)
    if(se == 0) se <- 1e-10 
    stat <- (p_hat - p0) / se
    
    if (alternative == "less") p_val <- pnorm(stat)
    else if (alternative == "greater") p_val <- 1 - pnorm(stat)
    else p_val <- 2 * (1 - pnorm(abs(stat)))
    
    return(list(statistic = stat, p.value = p_val, method = "Z-test Approx (Large n)"))
  } else {
    res <- binom.test(x_count, n, p = p0, alternative = alternative)
    return(list(statistic = res$statistic, p.value = res$p.value, method = "Exact Binomial (Small n)"))
  }
}

#' @param x1 Vector 1 o lista summary.
#' @param x2 Vector 2 o lista summary.
#' @param sigma1 Sigma conocido pob 1 (opcional).
#' @param sigma2 Sigma conocido pob 2 (opcional).
#' @param var_equal Asumir varianzas iguales (default FALSE).
#' @param alternative "two.sided", "less", "greater".
#' @export
test_diff_means <- function(x1, x2, sigma1 = NULL, sigma2 = NULL, var_equal = FALSE, alternative = "two.sided") {
  # Handle input x1
  if (is.numeric(x1) && length(x1) > 1) {
    n1 <- length(x1); m1 <- mean(x1); s1 <- sd(x1)
  } else {
    n1 <- x1$n; m1 <- x1$mean; s1 <- x1$sd
  }
  
  # Handle input x2
  if (is.numeric(x2) && length(x2) > 1) {
    n2 <- length(x2); m2 <- mean(x2); s2 <- sd(x2)
  } else {
    n2 <- x2$n; m2 <- x2$mean; s2 <- x2$sd
  }
  
  diff_mean <- m1 - m2
  
  if (!is.null(sigma1) && !is.null(sigma2)) {
    # Known Variances
    se <- sqrt(sigma1^2/n1 + sigma2^2/n2)
    stat <- diff_mean / se
    dist <- "norm"
    method <- "Z-test (Known Variances)"
  } else {
    if (n1 + n2 > 30) {
      # Large Samples
      se <- sqrt(s1^2/n1 + s2^2/n2)
      stat <- diff_mean / se
      dist <- "norm"
      method <- "Z-test Approx (Large Samples)"
    } else {
      # Small Samples
      if (var_equal) {
        df <- n1 + n2 - 2
        sp_sq <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / df
        se <- sqrt(sp_sq * (1/n1 + 1/n2))
        stat <- diff_mean / se
        dist <- "t"
        method <- "T-test (Small n, Pooled)"
      } else {
        # Welch
        num <- (s1^2/n1 + s2^2/n2)^2
        den <- ((s1^2/n1)^2 / (n1 - 1)) + ((s2^2/n2)^2 / (n2 - 1))
        df <- num / den
        se <- sqrt(s1^2/n1 + s2^2/n2)
        stat <- diff_mean / se
        dist <- "t"
        method <- "T-test (Small n, Welch)"
      }
    }
  }
  
  if (dist == "norm") {
    if (alternative == "less") p_val <- pnorm(stat)
    else if (alternative == "greater") p_val <- 1 - pnorm(stat)
    else p_val <- 2 * (1 - pnorm(abs(stat)))
  } else {
    if (alternative == "less") p_val <- pt(stat, df = df)
    else if (alternative == "greater") p_val <- 1 - pt(stat, df = df)
    else p_val <- 2 * (1 - pt(abs(stat), df = df))
  }
  
  return(list(statistic = stat, p.value = p_val, method = method, alternative = alternative))
}

#' @param x1 Éxitos muestra 1.
#' @param n1 Tamaño muestra 1.
#' @param x2 Éxitos muestra 2.
#' @param n2 Tamaño muestra 2.
#' @param alternative "two.sided", "less", "greater".
#' @export
test_diff_proportions <- function(x1, n1, x2, n2, alternative = "two.sided") {
  p1_hat <- x1 / n1
  p2_hat <- x2 / n2
  
  # Using separate variances for SE calculation
  se <- sqrt( (p1_hat * (1 - p1_hat) / n1) + (p2_hat * (1 - p2_hat) / n2) )
  
  if (se == 0) se <- 1e-10
  stat <- (p1_hat - p2_hat) / se
  
  if (alternative == "less") p_val <- pnorm(stat)
  else if (alternative == "greater") p_val <- 1 - pnorm(stat)
  else p_val <- 2 * (1 - pnorm(abs(stat)))
  
  return(list(statistic = stat, p.value = p_val, method = "Z-test (Diff Proportions)"))
}

#' @param observed Vector de frecuencias observadas.
#' @param expected Probabilidades esperadas O frecuencias esperadas.
#' @param num_params Número de parámetros estimados (m). Default 0.
#' @param alpha Nivel de significancia. Default 0.05.
#' @export
test_goodness_fit <- function(observed, expected, num_params = 0, alpha = 0.05) {
  
  # Setup data
  n <- sum(observed)
  k_initial <- length(observed)
  
  if (abs(sum(expected) - 1) < 0.01) {
    expected_freqs <- expected * n
  } else {
    expected_freqs <- expected
    if (abs(sum(expected_freqs) - n) > 0.1) {
      warning("Expected frequencies do not sum to N. Re-normalizing...")
      expected_freqs <- (expected_freqs / sum(expected_freqs)) * n
    }
  }
  
  if (is.null(names(observed))) {
    bin_names <- as.character(0:(k_initial-1))
    if (length(bin_names) != k_initial) bin_names <- paste0("Bin", 1:k_initial)
  } else {
    bin_names <- names(observed)
  }
  
  # Dataframe for bin merging logic
  df_work <- data.frame(
    Name = bin_names,
    O = observed,
    E = expected_freqs,
    stringsAsFactors = FALSE
  )
  
  # Merge bins where Expected < 5
  merged_flag <- TRUE
  while (merged_flag) {
    merged_flag <- FALSE
    rows <- nrow(df_work)
    
    if (rows <= 1) break 
    
    # Check Tail: Last bin
    if (df_work$E[rows] < 5) {
      df_work$Name[rows-1] <- paste0(df_work$Name[rows-1], "+", df_work$Name[rows])
      df_work$O[rows-1]    <- df_work$O[rows-1] + df_work$O[rows]
      df_work$E[rows-1]    <- df_work$E[rows-1] + df_work$E[rows]
      df_work <- df_work[-rows, ]
      merged_flag <- TRUE
      next
    }
    
    # Check Head: First bin
    if (df_work$E[1] < 5) {
      df_work$Name[2] <- paste0(df_work$Name[1], "+", df_work$Name[2])
      df_work$O[2]    <- df_work$O[1] + df_work$O[2]
      df_work$E[2]    <- df_work$E[1] + df_work$E[2]
      df_work <- df_work[-1, ]
      merged_flag <- TRUE
      next
    }
    
    # Check Middle bins
    min_E <- min(df_work$E)
    if (min_E < 5) {
      idx <- which.min(df_work$E)
      # Determine merger target
      if (idx == 1) target <- 2
      else if (idx == nrow(df_work)) target <- idx - 1
      else {
        if (df_work$E[idx-1] < df_work$E[idx+1]) target <- idx - 1 else target <- idx + 1
      }
      
      dest <- min(idx, target)
      src  <- max(idx, target)
      
      df_work$Name[dest] <- paste0(df_work$Name[dest], "&", df_work$Name[src])
      df_work$O[dest]    <- df_work$O[dest] + df_work$O[src]
      df_work$E[dest]    <- df_work$E[dest] + df_work$E[src]
      df_work <- df_work[-src, ]
      merged_flag <- TRUE
    }
  }
  
  # Calc Stats
  df_work$Chi_Term <- (df_work$O - df_work$E)^2 / df_work$E
  chi_sq_calc <- sum(df_work$Chi_Term)
  
  k_final <- nrow(df_work)
  df <- k_final - num_params - 1
  
  if (df <= 0) {
    warning("Degrees of freedom <= 0.")
    crit_val <- NA
    p_val <- NA
    conclusion <- "Error: Not enough bins."
  } else {
    crit_val <- qchisq(1 - alpha, df)
    p_val <- pchisq(chi_sq_calc, df, lower.tail = FALSE)
    
    if (chi_sq_calc < crit_val) {
      conclusion <- "Do NOT Reject H0 (Fit is Good)"
    } else {
      conclusion <- "Reject H0 (Fit is NOT Good)"
    }
  }
  
  # Print Report
  cat(" Calculation Table:\n")
  print(format(df_work, digits = 4), row.names = FALSE)
  cat("\n")
  
  cat(" Test Statistics:\n")
  cat("   Sum of (O-E)^2/E (Chi^2 Calc): ", round(chi_sq_calc, 4), "\n")
  cat("   Estimated Parameters (m):       ", num_params, "\n")
  cat("   Degrees of Freedom (k - m - 1):", df, "\n")
  cat("   Significance Level (alpha):     ", alpha, "\n")
  cat("   Critical Value (Chi^2 table):   ", round(crit_val, 4), "\n\n")
  
  cat(" Conclusion:\n")
  cat("   ", chi_sq_calc, ifelse(chi_sq_calc < crit_val, "<", ">"), crit_val, "\n")
  cat("   Result:", conclusion, "\n")
  cat("=======================================================\n")
  
  invisible(list(
    statistic = chi_sq_calc,
    parameter = c(df = df),
    p.value = p_val,
    method = "Chi-square Goodness of Fit (Explicit)",
    data.name = deparse(substitute(observed)),
    table = df_work
  ))
}

#' @param table_matrix Matriz o tabla de contingencia.
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

#' @param x1 Vector numérico muestra 1.
#' @param x2 Vector numérico muestra 2.
#' @param alternative "two.sided", "less", "greater".
#' @export
test_var_ratio <- function(x1, x2, alternative = "two.sided") {
  res <- var.test(x1, x2, alternative = alternative)
  return(list(statistic = res$statistic, df = res$parameter, p.value = res$p.value, conf.int = res$conf.int))
}