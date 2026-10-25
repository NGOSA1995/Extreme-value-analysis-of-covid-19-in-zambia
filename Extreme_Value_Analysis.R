# ==========================================================
# GLOBAL OPTIONS & LIBRARIES
# ==========================================================
options(device.ask.default = FALSE)

library(tidyverse)
library(lubridate)
library(extRemes)
library(tseries)
library(evd)
library(boot)
library(xtable)
library(scales)
library(goftest) 
library(ggplot2)
library(knitr)
library(dplyr)
library(purrr)

# Global vector to store raw LaTeX table strings for file export
latex_tables_export <- character()

# ==========================================================
# 1. PUBLICATION-READY LATEX TABLE GENERATOR
# ==========================================================
print_latex_table <- function(df, caption, label, digits = 2) {
  
  df_clean <- df
  
  # Clean column headers for LaTeX rendering
  colnames(df_clean) <- gsub("_", " ", colnames(df_clean))
  colnames(df_clean) <- gsub("Theta", "$\\\\theta$", colnames(df_clean))
  colnames(df_clean) <- gsub("Lambda", "$\\\\lambda$", colnames(df_clean))
  colnames(df_clean) <- gsub("Sigma", "$\\\\sigma$", colnames(df_clean))
  colnames(df_clean) <- gsub("Xi", "$\\\\xi$", colnames(df_clean))
  colnames(df_clean) <- gsub("Mu", "$\\\\mu$", colnames(df_clean))
  
  # Format Parameter column if present
  if ("Parameter" %in% colnames(df_clean)) {
    df_clean$Parameter <- gsub("\\\\mu", "$\\\\mu$", df_clean$Parameter)
    df_clean$Parameter <- gsub("\\\\sigma", "$\\\\sigma$", df_clean$Parameter)
    df_clean$Parameter <- gsub("\\\\xi", "$\\\\xi$", df_clean$Parameter)
    df_clean$Parameter <- gsub("\\\\theta", "$\\\\theta$", df_clean$Parameter)
  }
  
  align_vec <- c("l", sapply(df_clean, function(col) if(is.numeric(col)) "r" else "l"))
  
  xtab <- xtable(
    df_clean, 
    caption = caption, 
    label   = label, 
    align   = align_vec,
    digits  = digits
  )
  
  xtab_out <- capture.output(
    print(
      xtab,
      include.rownames       = FALSE,
      booktabs               = TRUE,
      comment                = FALSE,
      caption.placement      = "top",
      sanitize.text.function = function(x) x
    )
  )
  
  clean_str <- paste(xtab_out, collapse = "\n")
  cat("\n", clean_str, "\n\n")
  
  latex_tables_export <<- c(latex_tables_export, clean_str, "\n\n\\clearpage\n\n")
}

# ==========================================================
# 2. HELPER FUNCTIONS: GEV QUANTILE & GOF METRICS
# ==========================================================
gev_quantile <- function(p, loc, scale, shape) {
  if (is.null(p) || length(p) == 0 || any(is.na(p))) {
    return(rep(NA_real_, length(p)))
  }
  res <- numeric(length(p))
  for (i in seq_along(p)) {
    pi <- p[i]
    if (is.na(pi) || pi <= 0 || pi >= 1) {
      res[i] <- NA_real_
    } else if (abs(shape) < 1e-8) {
      res[i] <- loc - scale * log(-log(pi))
    } else {
      res[i] <- loc + (scale / shape) * ((-log(pi))^(-shape) - 1)
    }
  }
  return(res)
}

ks_statistic <- function(F_values) {
  F_values <- sort(F_values)
  n <- length(F_values)
  i <- seq_len(n)
  D_plus  <- max(i / n - F_values)
  D_minus <- max(F_values - (i - 1) / n)
  max(D_plus, D_minus)
}

ad_statistic <- function(F_values) {
  F_values <- sort(F_values)
  n <- length(F_values)
  F_values <- pmin(pmax(F_values, 1e-10), 1 - 1e-10)
  i <- seq_len(n)
  AD <- -n - mean((2 * i - 1) * (log(F_values) + log(1 - rev(F_values))))
  return(AD)
}

gof_boot_gev <- function(x_pos, B = 2000) {
  N <- length(x_pos)
  if (N < 5) return(c(KS_p = NA_real_, AD_p = NA_real_))
  
  fit_obs <- tryCatch(fevd(x_pos, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit_obs)) return(c(KS_p = NA_real_, AD_p = NA_real_))
  
  pars_obs <- distill(fit_obs)
  L_obs <- as.numeric(pars_obs["location"])
  S_obs <- as.numeric(pars_obs["scale"])
  G_obs <- as.numeric(pars_obs["shape"])
  
  if (is.na(S_obs) || S_obs <= 0) return(c(KS_p = NA_real_, AD_p = NA_real_))
  
  F_obs <- pevd(x_pos, loc = L_obs, scale = S_obs, shape = G_obs, type = "GEV")
  F_obs <- pmin(pmax(F_obs, 1e-10), 1 - 1e-10)
  
  ks_obs <- ks_statistic(F_obs)
  ad_obs <- ad_statistic(F_obs)
  
  ks_boot <- numeric(B)
  ad_boot <- numeric(B)
  valid_count <- 0
  
  for (i in seq_len(B)) {
    sim_data <- tryCatch(evd::rgev(N, loc = L_obs, scale = S_obs, shape = G_obs), error = function(e) NULL)
    if (is.null(sim_data) || any(!is.finite(sim_data))) next
    
    fit_sim <- tryCatch(fevd(sim_data, type = "GEV", method = "Lmoments"), error = function(e) NULL)
    if (is.null(fit_sim)) next
    
    p_sim <- distill(fit_sim)
    L_sim <- as.numeric(p_sim["location"])
    S_sim <- as.numeric(p_sim["scale"])
    G_sim <- as.numeric(p_sim["shape"])
    
    if (is.na(S_sim) || S_sim <= 0) next
    
    F_sim <- pevd(sim_data, loc = L_sim, scale = S_sim, shape = G_sim, type = "GEV")
    F_sim <- pmin(pmax(F_sim, 1e-10), 1 - 1e-10)
    
    ks_star <- tryCatch(ks_statistic(F_sim), error = function(e) NA_real_)
    ad_star <- tryCatch(ad_statistic(F_sim), error = function(e) NA_real_)
    
    if (is.finite(ks_star) && is.finite(ad_star)) {
      valid_count <- valid_count + 1
      ks_boot[valid_count] <- ks_star
      ad_boot[valid_count] <- ad_star
    }
  }
  
  if (valid_count < 100) return(c(KS_p = NA_real_, AD_p = NA_real_))
  
  ks_boot <- ks_boot[seq_len(valid_count)]
  ad_boot <- ad_boot[seq_len(valid_count)]
  
  ks_pval <- (sum(ks_boot >= ks_obs) + 1) / (valid_count + 1)
  ad_pval <- (sum(ad_boot >= ad_obs) + 1) / (valid_count + 1)
  
  return(c(KS_p = ks_pval, AD_p = ad_pval))
}

zero_adjusted_return_level <- function(T_years, loc, scale, shape, p_zero, blocks_per_year) {
  p_unconditional <- 1 - 1 / (blocks_per_year * T_years)
  p_positive      <- 1 - p_zero
  p_gev <- (p_unconditional - p_zero) / p_positive
  p_gev <- max(min(p_gev, 1 - 1e-10), 1e-10)
  gev_quantile(p = p_gev, loc = loc, scale = scale, shape = shape)
}

# ==========================================================
# 3. DATA PREPARATION & STRUCTURE
# ==========================================================
zambia_data <- Cov_Za %>%
  mutate(
    Day = as.Date(Day),
    cases = as.numeric(Daily.new.confirmed.cases.of.COVID.19.per.million.people)
  ) %>%
  filter(
    Entity == "Zambia",
    Day >= as.Date("2020-01-04"),
    Day <= as.Date("2023-12-31"),
    !is.na(cases)
  ) %>%
  arrange(Day)

cases_vector <- zambia_data$cases

complete_cases <- zambia_data %>%
  filter(!is.na(cases)) %>%
  mutate(Day = as.Date(Day)) %>%
  arrange(Day)

active_cases <- complete_cases %>%
  filter(cases > 0) %>%
  pull(cases)

n_total_obs    <- nrow(complete_cases)
n_zero_obs     <- sum(complete_cases$cases == 0, na.rm = TRUE)
n_positive_obs <- length(active_cases)

study_days  <- as.numeric(max(complete_cases$Day) - min(complete_cases$Day)) + 1
study_years <- study_days / 365.25

# ==========================================================
# 4. BLOCK MAXIMA EVALUATION (TABLE 1)
# ==========================================================
results_list <- list()
block_sizes  <- c(7, 14, 21, 28)

final_maxima_all <- NULL
final_maxima     <- NULL
final_fit        <- NULL
final_p_zero     <- NULL
final_p_positive <- NULL
final_n_blocks   <- NULL

set.seed(123)

for (b in block_sizes) {
  n_full_blocks <- floor(length(cases_vector) / b)
  if (n_full_blocks < 5) next
  
  trimmed_data <- cases_vector[1:(n_full_blocks * b)]
  data_matrix  <- matrix(trimmed_data, nrow = b, byrow = FALSE)
  
  block_maxima_all <- apply(data_matrix, 2, max, na.rm = TRUE)
  
  total_blocks    <- length(block_maxima_all)
  zero_blocks     <- sum(block_maxima_all == 0)
  positive_blocks <- sum(block_maxima_all > 0)
  
  p_zero     <- zero_blocks / total_blocks
  p_positive <- positive_blocks / total_blocks
  
  positive_maxima <- block_maxima_all[block_maxima_all > 0]
  if (length(positive_maxima) < 5) next
  
  tryCatch({
    gev_fit <- fevd(positive_maxima, type = "GEV", method = "Lmoments")
    params  <- distill(gev_fit)
    
    L <- as.numeric(params["location"])
    S <- as.numeric(params["scale"])
    G <- as.numeric(params["shape"])
    
    N        <- length(positive_maxima)
    p_seq    <- ppoints(N)
    q_theory <- qevd(p_seq, loc = L, scale = S, shape = G, type = "GEV")
    rmse_positive <- sqrt(mean((sort(positive_maxima) - q_theory)^2))
    
    gof_pvals <- gof_boot_gev(positive_maxima, B = 2000)
    
    diffs     <- diff(block_maxima_all)
    signs     <- sign(diffs)[sign(diffs) != 0]
    runs_pval <- if (length(signs) > 5) {
      tryCatch(tseries::runs.test(as.factor(signs > 0))$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    
    results_list[[as.character(b)]] <- data.frame(
      `Block Size` = b,
      `n Blks`     = total_blocks,
      `Zero Blks`  = zero_blocks,
      `Pos Blks`   = positive_blocks,
      `P(Zero)`    = p_zero,
      `P(Pos)`     = p_positive,
      `RMSE+`      = rmse_positive,
      `Runs (p)`   = runs_pval,
      `KS (pboot)` = as.numeric(gof_pvals["KS_p"]),
      `AD (pboot)` = as.numeric(gof_pvals["AD_p"]),
      `mu`         = L,
      `sigma`      = S,
      `xi`         = G,
      check.names  = FALSE
    )
    
    if (b == 28) {
      final_maxima_all <<- block_maxima_all
      final_maxima     <<- positive_maxima
      final_fit        <<- gev_fit
      final_p_zero     <<- p_zero
      final_p_positive <<- p_positive
      final_n_blocks   <<- total_blocks
    }
  }, error = function(e) {
    message(paste("Error at block size", b, ":", e$message))
  })
}

table1_results <- do.call(rbind, results_list)
rownames(table1_results) <- NULL

print_latex_table(
  table1_results,
  caption = "Block maxima model evaluation retaining zero-valued blocks.",
  label   = "tab:table1_block_size_eval",
  digits  = c(0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3)
)

# ==========================================================
# 5. GEV L-MOMENTS PARAMETERS (TABLE 4)
# ==========================================================
pars_lmom <- distill(final_fit)
loc_hat   <- as.numeric(pars_lmom["location"])
scale_hat <- as.numeric(pars_lmom["scale"])
shape_hat <- as.numeric(pars_lmom["shape"])

boot_gev_pars <- function(data, indices) {
  d <- data[indices]
  fit <- tryCatch(fevd(d, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit)) return(c(location = NA, scale = NA, shape = NA))
  p <- distill(fit)
  c(location = as.numeric(p["location"]), scale = as.numeric(p["scale"]), shape = as.numeric(p["shape"]))
}

set.seed(123)
boot_pars_res <- boot(data = final_maxima, statistic = boot_gev_pars, R = 2000)

# FIX: Explicitly set matrix column names for boot result array
colnames(boot_pars_res$t) <- c("location", "scale", "shape")

table4_lmom_param <- data.frame(
  Parameter    = c("Location (\\mu)", "Scale (\\sigma)", "Shape (\\xi)"),
  Estimate     = c(loc_hat, scale_hat, shape_hat),
  `Std Error`  = apply(boot_pars_res$t, 2, sd, na.rm = TRUE),
  `Lower 95% CI` = apply(boot_pars_res$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  `Upper 95% CI` = apply(boot_pars_res$t, 2, quantile, probs = 0.975, na.rm = TRUE),
  check.names  = FALSE
)

print_latex_table(
  table4_lmom_param,
  caption = "GEV parameter estimates obtained via L-moments with bootstrap standard errors and 95\\% confidence intervals.",
  label   = "tab:table4_gev_lmom_uncertainty",
  digits  = 4
)

# ==========================================================
# 6. ZERO-ADJUSTED STATIONARY RETURN LEVELS (TABLE 5)
# ==========================================================
blocks_per_year      <- 365.25 / 28
return_periods_years <- c(2, 5, 10, 20)

point_rl <- sapply(return_periods_years, function(T) {
  zero_adjusted_return_level(T_years = T, loc = loc_hat, scale = scale_hat, shape = shape_hat, p_zero = final_p_zero, blocks_per_year = blocks_per_year)
})

boot_zero_adjusted_rl <- function(data, indices, periods, blocks_per_year) {
  d <- data[indices]
  p_zero_b     <- mean(d == 0)
  d_positive   <- d[d > 0]
  if (length(d_positive) < 5 || (1 - p_zero_b) <= 0) return(rep(NA_real_, length(periods)))
  fit_b <- tryCatch(fevd(d_positive, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit_b)) return(rep(NA_real_, length(periods)))
  pars_b <- tryCatch(distill(fit_b), error = function(e) NULL)
  if (is.null(pars_b)) return(rep(NA_real_, length(periods)))
  s_b <- as.numeric(pars_b["scale"])
  if (is.na(s_b) || s_b <= 0) return(rep(NA_real_, length(periods)))
  
  sapply(periods, function(T) {
    zero_adjusted_return_level(T_years = T, loc = as.numeric(pars_b["location"]), scale = s_b, shape = as.numeric(pars_b["shape"]), p_zero = p_zero_b, blocks_per_year = blocks_per_year)
  })
}

set.seed(123)
boot_rl_res <- boot(data = final_maxima_all, statistic = boot_zero_adjusted_rl, R = 2000, periods = return_periods_years, blocks_per_year = blocks_per_year)

table5_return_levels <- data.frame(
  `Return Period Years` = paste0(return_periods_years, "-year"),
  Estimate               = point_rl,
  `Lower 95% CI`         = apply(boot_rl_res$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  `Upper 95% CI`         = apply(boot_rl_res$t, 2, quantile, probs = 0.975, na.rm = TRUE),
  check.names            = FALSE
)

print_latex_table(
  table5_return_levels,
  caption = "Zero-adjusted stationary GEV return level estimates for Zambia based on 28-day block maxima.",
  label   = "tab:table5_zero_adjusted_return_levels",
  digits  = 2
)

# ==========================================================
# 7. REVISED SECTION 7: QUADRATIC GEV SENSITIVITY & TABLE 6
# ==========================================================

n_blocks_pos <- length(final_maxima)
time_cov_raw <- 1:n_blocks_pos
time_cov     <- as.numeric(scale(time_cov_raw))
time_cov_sq  <- time_cov^2
data_cov_pos <- data.frame(time_cov = time_cov, time_cov_sq = time_cov_sq)

fit_stat_mle <- fevd(final_maxima, type = "GEV", method = "MLE")
fit_lin_mle  <- fevd(final_maxima, type = "GEV", location.fun = ~ time_cov, data = data_cov_pos, method = "MLE")
fit_quad_mle <- fevd(final_maxima, type = "GEV", location.fun = ~ time_cov + time_cov_sq, data = data_cov_pos, method = "MLE")

models_mle  <- list(fit_stat_mle, fit_lin_mle, fit_quad_mle)
model_names <- c("Stationary", "Linear Trend", "Quadratic Trend")

loglik_values <- sapply(models_mle, function(x) -as.numeric(x$results$value))
aic_values    <- sapply(models_mle, function(x) as.numeric(summary(x, silent = TRUE)$AIC))
bic_values    <- sapply(models_mle, function(x) as.numeric(summary(x, silent = TRUE)$BIC))

n_parameters <- c(3, 4, 5)
aicc_values  <- aic_values + (2 * n_parameters * (n_parameters + 1)) / (n_blocks_pos - n_parameters - 1)

table2_model_comp <- data.frame(
  Model      = model_names,
  `Log Lik`  = loglik_values,
  AIC        = aic_values,
  AICc       = aicc_values,
  BIC        = bic_values,
  Parameters = n_parameters,
  check.names = FALSE
)

print_latex_table(
  table2_model_comp,
  caption = "Comparison of stationary and non-stationary Generalised Extreme Value models fitted to 28-day block maxima.",
  label   = "tab:table2_gev_model_comp",
  digits  = c(0, 0, 2, 2, 2, 2, 0)
)

lr_stat_lin  <- 2 * (loglik_values[2] - loglik_values[1])
lr_lin_quad  <- 2 * (loglik_values[3] - loglik_values[2])
lr_stat_quad <- 2 * (loglik_values[3] - loglik_values[1])

table3_lrt <- data.frame(
  Comparison = c("Stationary vs Linear Trend", "Linear Trend vs Quadratic Trend", "Stationary vs Quadratic Trend"),
  LR         = c(lr_stat_lin, lr_lin_quad, lr_stat_quad),
  df         = c(1, 1, 2),
  `p value`  = c(
    1 - pchisq(lr_stat_lin, df = 1),
    1 - pchisq(lr_lin_quad, df = 1),
    1 - pchisq(lr_stat_quad, df = 2)
  ),
  check.names = FALSE
)

print_latex_table(
  table3_lrt,
  caption = "Likelihood-ratio tests for nested Generalised Extreme Value models fitted to 28-day block maxima.",
  label   = "tab:table3_lrt",
  digits  = c(0, 0, 3, 0, 3)
)

# 1. Primary L-moment Stationary Return Levels (from Table 5)
pars_lmom <- distill(final_fit)
loc_hat   <- as.numeric(pars_lmom["location"])
scale_hat <- as.numeric(pars_lmom["scale"])
shape_hat <- as.numeric(pars_lmom["shape"])

return_periods_years <- c(2, 5, 10, 20)
blocks_per_year      <- 365.25 / 28

point_rl_lmom <- sapply(return_periods_years, function(T_yr) {
  zero_adjusted_return_level(
    T_years         = T_yr,
    loc             = loc_hat,
    scale           = scale_hat,
    shape           = shape_hat,
    p_zero          = final_p_zero,
    blocks_per_year = blocks_per_year
  )
})

# 2. Extract Quadratic MLE Parameter Estimates
required_names <- c("mu0", "mu1", "mu2", "scale", "shape")
quad_pars      <- fit_quad_mle$results$par[required_names]

beta0    <- unname(quad_pars["mu0"])
beta1    <- unname(quad_pars["mu1"])
beta2    <- unname(quad_pars["mu2"])
sig_quad <- unname(quad_pars["scale"])
xi_quad  <- unname(quad_pars["shape"])

cov_mat <- tryCatch(solve(fit_quad_mle$results$hessian), error = function(e) NULL)
se_quad <- if (!is.null(cov_mat) && all(diag(cov_mat) > 0)) sqrt(diag(cov_mat))[required_names] else rep(NA_real_, 5)

quad_parameter_table <- data.frame(
  Parameter = c(
    "Location Intercept (\\mu_0)",
    "Linear Location Coeff. (\\mu_1)",
    "Quadratic Location Coeff. (\\mu_2)",
    "Scale (\\sigma)",
    "Shape (\\xi)"
  ),
  Estimate     = as.numeric(quad_pars),
  `Std Error`  = as.numeric(se_quad),
  `Lower 95% CI` = as.numeric(quad_pars - 1.96 * se_quad),
  `Upper 95% CI` = as.numeric(quad_pars + 1.96 * se_quad),
  check.names  = FALSE
)

# 3. Compute End-of-Period Quadratic Return Levels
mu_quad     <- beta0 + beta1 * time_cov + beta2 * time_cov_sq
mu_quad_end <- tail(mu_quad, 1)

rl_quad_end <- sapply(return_periods_years, function(T_yr) {
  zero_adjusted_return_level(
    T_years         = T_yr,
    loc             = mu_quad_end,
    scale           = sig_quad,
    shape           = xi_quad,
    p_zero          = final_p_zero,
    blocks_per_year = blocks_per_year
  )
})

# 4. Parametric Bootstrap for End-of-Period Quadratic CIs
set.seed(123)
B_boot <- 2000

boot_quad_rl <- matrix(
  NA_real_, 
  nrow = B_boot, 
  ncol = length(return_periods_years),
  dimnames = list(NULL, paste0(return_periods_years, "-year"))
)

for (b in seq_len(B_boot)) {
  simulated_maxima <- tryCatch(
    evd::rgev(n = n_blocks_pos, loc = mu_quad, scale = sig_quad, shape = xi_quad),
    error = function(e) NULL
  )
  if (is.null(simulated_maxima) || any(!is.finite(simulated_maxima))) next
  
  fit_b <- tryCatch(
    suppressWarnings(
      extRemes::fevd(
        x = simulated_maxima, type = "GEV",
        location.fun = ~ time_cov + time_cov_sq,
        data = data_cov_pos, method = "MLE"
      )
    ),
    error = function(e) NULL
  )
  if (is.null(fit_b)) next
  
  if (!identical(as.integer(fit_b$results$convergence), 0L)) next
  
  pb <- fit_b$results$par
  if (is.null(pb) || !all(required_names %in% names(pb))) next
  
  pb <- pb[required_names]
  if (any(!is.finite(pb)) || pb["scale"] <= 0) next
  
  mu_b <- pb["mu0"] + pb["mu1"] * time_cov + pb["mu2"] * time_cov_sq
  
  rl_b <- sapply(return_periods_years, function(T_yr) {
    zero_adjusted_return_level(
      T_years         = T_yr,
      loc             = tail(mu_b, 1),
      scale           = unname(pb["scale"]),
      shape           = unname(pb["shape"]),
      p_zero          = final_p_zero,
      blocks_per_year = blocks_per_year
    )
  })
  
  if (all(is.finite(rl_b))) {
    boot_quad_rl[b, ] <- as.numeric(rl_b)
  }
}

valid_rl <- apply(boot_quad_rl, 1, function(x) all(is.finite(x)))
boot_quad_rl_valid <- boot_quad_rl[valid_rl, , drop = FALSE]

quad_rl_lower <- apply(boot_quad_rl_valid, 2, quantile, probs = 0.025, names = FALSE)
quad_rl_upper <- apply(boot_quad_rl_valid, 2, quantile, probs = 0.975, names = FALSE)

# 5. Construct Synchronized Table 6 Data Frame
table6_sensitivity <- data.frame(
  `Return Period`                = paste0(return_periods_years, "-year"),
  `Stationary L-moment Estimate` = as.numeric(point_rl_lmom),
  `Quadratic End Estimate`       = as.numeric(rl_quad_end),
  `Quadratic Lower 95% CI`       = as.numeric(quad_rl_lower),
  `Quadratic Upper 95% CI`       = as.numeric(quad_rl_upper),
  check.names                    = FALSE
)

# ==========================================================
# EXPLICIT PRINTING & EXPORTING OF TABLES
# ==========================================================

# A. Print raw R Data Frame to Console
cat("\n=== DATA FRAME OUTPUT (TABLE 6) ===\n")
print(table6_sensitivity)

# B. Print Formatted LaTeX Code to Console and Store in latex_tables_export
cat("\n=== LATEX TABLE OUTPUT (TABLE 6) ===\n")
print_latex_table(
  table6_sensitivity,
  caption = paste(
    "Stationary L-moment and quadratic MLE zero-adjusted GEV",
    "return levels (cases per million)."
  ),
  label  = "tab:table6_sensitivity",
  digits = 2
)

# C. Export CSV directly to disk
write.csv(table6_sensitivity, "Table6_Quadratic_Return_Level_Sensitivity.csv", row.names = FALSE)

# ==========================================================
# FORCE PRINTING OF TABLE 4 (QUADRATIC PARAMETERS) & TABLE 6
# ==========================================================

# 1. DISPLAY TABLE 4 IN CONSOLE
cat("\n==========================================================\n")
cat("TABLE 4: QUADRATIC GEV PARAMETER ESTIMATES\n")
cat("==========================================================\n")
print(quad_parameter_table) # Prints R Data Frame
print_latex_table(           # Prints formatted LaTeX
  quad_parameter_table,
  caption = "Quadratic GEV maximum-likelihood parameter estimates.",
  label   = "tab:quadratic_gev_parameters",
  digits  = 4
)

# 2. DISPLAY TABLE 6 IN CONSOLE
cat("\n==========================================================\n")
cat("TABLE 6: RETURN LEVEL SENSITIVITY ANALYSIS\n")
cat("==========================================================\n")
print(table6_sensitivity)  # Prints R Data Frame
print_latex_table(           # Prints formatted LaTeX
  table6_sensitivity,
  caption = "Stationary L-moment and quadratic MLE return levels.",
  label   = "tab:table6_sensitivity",
  digits  = 2
)

# ==========================================================
# RETURN LEVEL CURVES & SENSITIVITY VISUALIZATION (ggplot2)
# ==========================================================

# Smooth Curve Grid
T_grid <- seq(1.05, 25, length.out = 300)

rl_stat_grid <- sapply(T_grid, function(T_yr) {
  zero_adjusted_return_level(
    T_years         = T_yr,
    loc             = loc_hat,
    scale           = scale_hat,
    shape           = shape_hat,
    p_zero          = final_p_zero,
    blocks_per_year = blocks_per_year
  )
})

boot_stat_grid <- matrix(NA_real_, nrow = 2000, ncol = length(T_grid))
for (b in seq_len(2000)) {
  # FIX: Use safe column indexing for the boot parameter matrix
  loc_b   <- boot_pars_res$t[b, "location"]
  scale_b <- boot_pars_res$t[b, "scale"]
  shape_b <- boot_pars_res$t[b, "shape"]
  
  if (is.finite(scale_b) && scale_b > 0) {
    boot_stat_grid[b, ] <- sapply(T_grid, function(T_yr) {
      zero_adjusted_return_level(
        T_years         = T_yr,
        loc             = loc_b,
        scale           = scale_b,
        shape           = shape_b,
        p_zero          = final_p_zero,
        blocks_per_year = blocks_per_year
      )
    })
  }
}

stat_lower_grid <- apply(boot_stat_grid, 2, quantile, probs = 0.025, na.rm = TRUE)
stat_upper_grid <- apply(boot_stat_grid, 2, quantile, probs = 0.975, na.rm = TRUE)

df_stat_curve <- data.frame(
  Return_Period = T_grid,
  Return_Level  = rl_stat_grid,
  Lower_CI      = stat_lower_grid,
  Upper_CI      = stat_upper_grid,
  Model         = "Stationary L-moment (Baseline)"
)

rl_quad_grid <- sapply(T_grid, function(T_yr) {
  zero_adjusted_return_level(
    T_years         = T_yr,
    loc             = mu_quad_end,
    scale           = sig_quad,
    shape           = xi_quad,
    p_zero          = final_p_zero,
    blocks_per_year = blocks_per_year
  )
})

df_quad_curve <- data.frame(
  Return_Period = T_grid,
  Return_Level  = rl_quad_grid,
  Lower_CI      = NA_real_,
  Upper_CI      = NA_real_,
  Model         = "Quadratic Trend (t = T_end)"
)

df_points <- data.frame(
  Return_Period = return_periods_years,
  Stationary    = as.numeric(point_rl_lmom),
  Quadratic     = as.numeric(rl_quad_end)
) %>%
  pivot_longer(cols = c(Stationary, Quadratic), names_to = "Model_Type", values_to = "Return_Level") %>%
  mutate(Model = if_else(Model_Type == "Stationary", "Stationary L-moment (Baseline)", "Quadratic Trend (t = T_end)"))

# FIGURE 2: Stationary Return Level Curve
p_fig2 <- ggplot(df_stat_curve, aes(x = Return_Period, y = Return_Level)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "#1F77B4", alpha = 0.20) +
  geom_line(color = "#1F77B4", linewidth = 1.1) +
  geom_point(data = filter(df_points, Model_Type == "Stationary"), 
             aes(x = Return_Period, y = Return_Level), 
             color = "#8C1D40", size = 3) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20), labels = c("1", "2", "5", "10", "20")) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Return Period (Years)",
    y = "Return Level (Cases per Million)",
    title = "Zero-Adjusted Stationary GEV Return Level Plot (28-day Block Maxima)",
    subtitle = "Solid curve represents point estimates; shaded region denotes 95% bootstrap CI."
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    axis.title       = element_text(face = "bold")
  )

print(p_fig2)
ggsave("Figure2_Stationary_GEV_Return_Levels.pdf", plot = p_fig2, width = 8, height = 5, dpi = 300)
ggsave("Figure2_Stationary_GEV_Return_Levels.png", plot = p_fig2, width = 8, height = 5, dpi = 300)

# FIGURE: Diagnostic Sensitivity Comparison
df_combined_curves <- bind_rows(df_stat_curve, df_quad_curve)

p_sensitivity <- ggplot() +
  geom_ribbon(data = df_stat_curve, aes(x = Return_Period, ymin = Lower_CI, ymax = Upper_CI), 
              fill = "#1F77B4", alpha = 0.15) +
  geom_line(data = df_combined_curves, aes(x = Return_Period, y = Return_Level, color = Model, linetype = Model), 
            linewidth = 1.1) +
  geom_point(data = df_points, aes(x = Return_Period, y = Return_Level, color = Model), size = 3) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20), labels = c("1", "2", "5", "10", "20")) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = c("Stationary L-moment (Baseline)" = "#1F77B4", "Quadratic Trend (t = T_end)" = "#D95F02")) +
  scale_linetype_manual(values = c("Stationary L-moment (Baseline)" = "solid", "Quadratic Trend (t = T_end)" = "dashed")) +
  labs(
    x = "Return Period (Years, Log Scale)",
    y = "Return Level (Cases per Million, Log Scale)",
    title = "Sensitivity of Return Level Extrapolation to Location Drift",
    subtitle = "Comparing Stationary L-moment baseline vs Non-Stationary Quadratic boundary fit",
    color = "Model Specification",
    linetype = "Model Specification"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    axis.title       = element_text(face = "bold")
  )

print(p_sensitivity)
ggsave("Figure_Sensitivity_Stationary_vs_Quadratic.pdf", plot = p_sensitivity, width = 8.5, height = 5.5, dpi = 300)
ggsave("Figure_Sensitivity_Stationary_vs_Quadratic.png", plot = p_sensitivity, width = 8.5, height = 5.5, dpi = 300)

# ==========================================================
# 8. POT ANALYSIS & THRESHOLD SENSITIVITY (TABLES 7, 8, AND 9)
# ==========================================================
percentiles <- seq(0.85, 0.95, by = 0.01)

u_complete <- quantile(complete_cases$cases, probs = percentiles, na.rm = TRUE, type = 7)
u_positive <- quantile(active_cases, probs = percentiles, na.rm = TRUE, type = 7)

evaluate_threshold <- function(u, percentile_label, data, n_total, years, run_length = 7) {
  exceedances <- data[data > u]
  raw_N        <- length(exceedances)
  if (raw_N == 0) return(NULL)
  
  declustered_series <- decluster(x = data, threshold = u, method = "runs", r = run_length)
  cluster_N <- sum(declustered_series > u, na.rm = TRUE)
  theta <- cluster_N / raw_N
  raw_rate_daily  <- raw_N / n_total
  lambda_c_annual <- cluster_N / years
  
  gpd_fit <- tryCatch({ fevd(x = declustered_series, threshold = u, type = "GP", method = "Lmoments") }, error = function(e) NULL)
  
  if (is.null(gpd_fit)) {
    sigma_hat <- NA_real_; xi_hat <- NA_real_
  } else {
    pars <- distill(gpd_fit)
    sigma_hat <- as.numeric(pars["scale"])
    xi_hat    <- as.numeric(pars["shape"])
  }
  
  data.frame(
    Percentile           = percentile_label,
    `Threshold (u)`      = u,
    `Raw N (N_u)`        = raw_N,
    `Clusters (N_c)`     = cluster_N,
    `Theta (\\theta)`    = theta,
    `Raw Prop (N_u/n_+)` = raw_rate_daily,
    `Cluster Rate (\\lambda_c)` = lambda_c_annual,
    `Scale (\\sigma)`    = sigma_hat,
    `Shape (\\xi)`       = xi_hat,
    check.names          = FALSE
  )
}

table7_positive <- map2_df(u_positive, paste0(round(percentiles * 100), "\\% positive"), ~ evaluate_threshold(.x, .y, active_cases, n_positive_obs, study_years, run_length = 7))
table7_complete <- map2_df(u_complete, paste0(round(percentiles * 100), "\\% complete"), ~ evaluate_threshold(.x, .y, complete_cases$cases, n_total_obs, study_years, run_length = 7))

table7_full <- bind_rows(table7_positive, table7_complete)

print_latex_table(
  table7_full,
  caption = "Threshold sensitivity analysis comparing positive-observation and complete-series percentiles.",
  label   = "tab:table7_threshold_sensitivity",
  digits  = c(0, 0, 2, 0, 0, 3, 3, 3, 3, 3)
)

u_selected <- 61.26
run_length <- 7

exceedance_data <- complete_cases %>%
  filter(cases > u_selected) %>%
  arrange(Day) %>%
  mutate(
    gap_days    = as.numeric(Day - lag(Day)),
    new_cluster = if_else(is.na(gap_days) | gap_days > run_length, 1, 0),
    cluster_id  = cumsum(new_cluster)
  )

table8_cluster_data <- exceedance_data %>%
  group_by(cluster_id) %>%
  summarise(
    `Cluster Start`   = as.character(min(Day)),
    `Cluster End`     = as.character(max(Day)),
    `Duration Days`   = as.numeric(max(Day) - min(Day)) + 1,
    `Cluster Size`    = n(),
    `Cluster Maximum` = max(cases),
    .groups           = "drop"
  )

print_latex_table(
  table8_cluster_data,
  caption = "Characteristics of independent exceedance clusters after runs declustering at the selected 90th percentile threshold.",
  label   = "tab:table8_cluster_characteristics",
  digits  = 2
)

run_lengths <- c(1, 2, 3, 4, 5, 7, 10, 14, 21)

table9_runs_sens <- map_df(run_lengths, function(r) {
  res <- evaluate_threshold(u_selected, paste0(r, " days"), active_cases, n_positive_obs, study_years, run_length = r)
  data.frame(
    `Run Length Days`      = r,
    `Raw Exceedances`      = res$`Raw N (N_u)`,
    `Independent Clusters` = res$`Clusters (N_c)`,
    Theta                  = res$`Theta (\\theta)`,
    `Annual Cluster Rate`  = res$`Cluster Rate (\\lambda_c)`,
    check.names          = FALSE
  )
})

print_latex_table(
  table9_runs_sens,
  caption = "Sensitivity of runs declustering to alternative run lengths at the selected 90th percentile threshold.",
  label   = "tab:table9_runs_sensitivity",
  digits  = c(0, 0, 0, 0, 2, 2)
)

# ==========================================================
# 9. GPD FITTING & POT RETURN LEVELS (TABLES 10 AND 11)
# ==========================================================
cluster_peaks <- table8_cluster_data$`Cluster Maximum`
cluster_sizes <- table8_cluster_data$`Cluster Size`

final_gpd_fit <- fevd(x = cluster_peaks, threshold = u_selected, type = "GP", method = "Lmoments")
gpd_pars      <- distill(final_gpd_fit)
sigma_gpd     <- as.numeric(gpd_pars["scale"])
xi_gpd        <- as.numeric(gpd_pars["shape"])
lambda_c      <- nrow(table8_cluster_data) / study_years

gpd_rl_estimates <- sapply(return_periods_years, function(T_yr) {
  if (abs(xi_gpd) < 1e-8) {
    u_selected + sigma_gpd * log(lambda_c * T_yr)
  } else {
    u_selected + (sigma_gpd / xi_gpd) * (((lambda_c * T_yr)^xi_gpd) - 1)
  }
})

boot_gpd_rl <- function(peaks, R = 5000) {
  res_mat <- matrix(NA_real_, nrow = R, ncol = length(return_periods_years))
  n_peaks <- length(peaks)
  for (i in 1:R) {
    samp <- sample(peaks, size = n_peaks, replace = TRUE)
    fit  <- tryCatch(fevd(samp, threshold = u_selected, type = "GP", method = "Lmoments"), error = function(e) NULL)
    if (is.null(fit)) next
    p    <- distill(fit)
    sig <- as.numeric(p["scale"])
    xi  <- as.numeric(p["shape"])
    if (is.na(sig) || sig <= 0) next
    
    res_mat[i, ] <- sapply(return_periods_years, function(T_yr) {
      if (abs(xi) < 1e-8) {
        u_selected + sig * log(lambda_c * T_yr)
      } else {
        u_selected + (sig / xi) * (((lambda_c * T_yr)^xi) - 1)
      }
    })
  }
  return(res_mat)
}

set.seed(123)
gpd_boot_rl_samples <- boot_gpd_rl(cluster_peaks, R = 5000)

table10_gpd_rl <- data.frame(
  `Return Period` = paste0(return_periods_years, "-year"),
  `Lower 95% CI`  = apply(gpd_boot_rl_samples, 2, quantile, probs = 0.025, na.rm = TRUE),
  Estimate        = gpd_rl_estimates,
  `Upper 95% CI`  = apply(gpd_boot_rl_samples, 2, quantile, probs = 0.975, na.rm = TRUE),
  check.names     = FALSE
)

print_latex_table(
  table10_gpd_rl,
  caption = "GPD return level estimates with 95\\% confidence intervals for selected return periods.",
  label   = "tab:table10_gpd_return_levels",
  digits  = 2
)

boot_pot_params <- function(peaks, sizes, R = 2000) {
  theta_boot <- numeric(R)
  sig_boot   <- numeric(R)
  xi_boot    <- numeric(R)
  n_clusters <- length(peaks)
  
  for (i in 1:R) {
    boot_indices <- sample(seq_along(sizes), size = n_clusters, replace = TRUE)
    samp_sizes   <- sizes[boot_indices]
    samp_peaks   <- peaks[boot_indices]
    
    theta_boot[i] <- length(samp_sizes) / sum(samp_sizes)
    
    fit <- tryCatch(fevd(samp_peaks, threshold = u_selected, type = "GP", method = "Lmoments"), error = function(e) NULL)
    if (!is.null(fit)) {
      p <- distill(fit)
      sig_boot[i] <- as.numeric(p["scale"])
      xi_boot[i]  <- as.numeric(p["shape"])
    }
  }
  
  data.frame(
    Parameter    = c("Extremal Index (\\theta)", "Scale (\\sigma)", "Shape (\\xi)"),
    Estimate     = c(length(peaks) / sum(sizes), sigma_gpd, xi_gpd),
    `Lower 95% CI` = c(quantile(theta_boot, 0.025, na.rm=TRUE), quantile(sig_boot, 0.025, na.rm=TRUE), quantile(xi_boot, 0.025, na.rm=TRUE)),
    `Upper 95% CI` = c(quantile(theta_boot, 0.975, na.rm=TRUE), quantile(sig_boot, 0.975, na.rm=TRUE), quantile(xi_boot, 0.975, na.rm=TRUE)),
    check.names  = FALSE
  )
}

set.seed(123)
table11_pot_uncertainty <- boot_pot_params(cluster_peaks, cluster_sizes, R = 2000)

print_latex_table(
  table11_pot_uncertainty,
  caption = "Bootstrap uncertainty estimates for the extremal index and Generalised Pareto Distribution parameters at the selected 90th percentile threshold.",
  label   = "tab:table11_pot_uncertainty",
  digits  = 2
)

# ==========================================================
# 10. SAVE CONSOLIDATED LATEX FILE & CSV OUTPUTS
# ==========================================================
writeLines(latex_tables_export, "All_Manuscript_Tables.tex")

write.csv(table1_results, "Table1_Block_Maxima_Evaluation.csv", row.names = FALSE)
write.csv(table2_model_comp, "Table2_GEV_Model_Comparison.csv", row.names = FALSE)
write.csv(table3_lrt, "Table3_Likelihood_Ratio_Tests.csv", row.names = FALSE)
write.csv(quad_parameter_table, "Table4_Quadratic_GEV_Parameters.csv", row.names = FALSE)
write.csv(table4_lmom_param, "Table4_GEV_Lmoments_Parameters.csv", row.names = FALSE)
write.csv(table5_return_levels, "Table5_GEV_Return_Levels.csv", row.names = FALSE)
write.csv(table6_sensitivity, "Table6_Quadratic_Return_Level_Sensitivity.csv", row.names = FALSE)
write.csv(table7_full, "Table7_Threshold_Sensitivity_Analysis.csv", row.names = FALSE)
write.csv(table8_cluster_data, "Table8_Cluster_Characteristics.csv", row.names = FALSE)
write.csv(table9_runs_sens, "Table9_Runs_Declustering_Sensitivity.csv", row.names = FALSE)
write.csv(table10_gpd_rl, "Table10_GPD_Return_Levels.csv", row.names = FALSE)
write.csv(table11_pot_uncertainty, "Table11_POT_Bootstrap_Uncertainty.csv", row.names = FALSE)

cat("\n==========================================================\n")
cat("SUCCESS: MANUSCRIPT TABLES & FIGURES SUCCESSFULLY PRINTED & EXPORTED\n")
cat("==========================================================\n")