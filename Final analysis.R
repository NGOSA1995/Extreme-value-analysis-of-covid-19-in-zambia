# ==========================================================
# 1. SETUP & DATA LOADING
# ==========================================================
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

# Data Preparation
# Assuming Cov_Za is your original dataframe
zambia_data <- Cov_Za %>%
  mutate(Day = as.Date(Day),
         cases = as.numeric(Daily.new.confirmed.cases.of.COVID.19.per.million.people)) %>%
  filter(Entity == "Zambia",
         Day >= "2020-01-04" & Day <= "2023-12-31",
         !is.na(cases)) %>%
  arrange(Day)

cases_vector <- zambia_data$cases

# ==========================================================
# 2. ANALYSIS LOOP (BLOCK MAXIMA SELECTION)
# ==========================================================
results <- data.frame()
block_sizes <- c(7, 14, 21, 28)

for (b in block_sizes) {
  n_full_blocks <- floor(length(cases_vector) / b)
  if (n_full_blocks < 5) next
  
  trimmed_data <- cases_vector[1:(n_full_blocks * b)]
  data_matrix <- matrix(trimmed_data, nrow = b, byrow = FALSE)
  block_maxima_temp <- apply(data_matrix, 2, max, na.rm = TRUE)
  block_maxima_temp <- block_maxima_temp[block_maxima_temp > 0]
  
  tryCatch({
    # Fit GEV model using L-moments
    gev_fit <- fevd(block_maxima_temp, type = "GEV", method = "Lmoments")
    params <- distill(gev_fit)
    
    L <- as.numeric(params["location"])
    S <- as.numeric(params["scale"])
    G <- as.numeric(params["shape"])
    
    # --- VALIDATION STATS ---
    p_seq <- ppoints(length(block_maxima_temp))
    q_theory <- qgev(p_seq, loc = L, scale = S, shape = G)
    
    # RMSE
    rmse_val <- sqrt(mean((sort(block_maxima_temp) - q_theory)^2))
    
    # Handle ties for the KS test by adding minuscule jitter
    jittered_maxima <- jitter(block_maxima_temp, amount = 1e-5)
    
    # KS Test
    ks_pval <- ks.test(jittered_maxima, "pgev", loc = L, scale = S, shape = G)$p.value
    
    # Anderson-Darling Test (via goftest package)
    ad_pval <- goftest::ad.test(block_maxima_temp, "pgev", loc = L, scale = S, shape = G)$p.value
    
    # --- 1. RUNS TEST FOR INDEPENDENCE (Directional Sign-Difference Method) ---
    diffs <- diff(block_maxima_temp)
    signs <- sign(diffs)[sign(diffs) != 0]
    runs_pval <- if(length(signs) > 5) {
      tseries::runs.test(as.factor(signs > 0))$p.value
    } else {
      NA
    }
    
    # --- 2. WILCOXON SIGNED-RANK TEST ---
    theo_median <- qgev(0.5, loc = L, scale = S, shape = G)
    wilcox_pval <- wilcox.test(block_maxima_temp, mu = theo_median)$p.value
    
    # Combine results into structured dataframe
    results <- rbind(results, data.frame(
      BlockSize = b, 
      nBlocks   = length(block_maxima_temp), 
      RMSE      = rmse_val,
      Runs_p    = runs_pval,
      Wilcox_p  = wilcox_pval,
      KS_p      = ks_pval,
      AD_p      = ad_pval,
      Loc       = L, 
      Scale     = S, 
      Shape     = G
    ))
    
    if(b == 28) {
      final_maxima <- block_maxima_temp
      final_fit <- gev_fit
      loc <- L; scale <- S; shape <- G
    }
  }, error = function(e) { message(paste("Error at block size", b, ":", e)) })
}

# ==========================================================
#  EXPORT TO LATEX
# ==========================================================
# Renaming columns cleanly to print nicely in LaTeX formulas (μ, σ, ξ)
colnames(results) <- c("Block Size", "n Blks", "RMSE", "Runs (p)", "Wilcox (p)", "KS (p)", "AD (p)", "mu", "sigma", "xi")

print(xtable(results, digits = 3), include.rownames = FALSE)
#=======================================================================================================
# Statistic function for non-parametric bootstrapping of GEV parameters via L-moments
boot_gev_lmom <- function(data, indices) {
  d <- data[indices]
  fit <- tryCatch(
    fevd(d, type = "GEV", method = "Lmoments"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(location = NA, scale = NA, shape = NA))
  p <- distill(fit)
  return(c(
    location = as.numeric(p["location"]),
    scale    = as.numeric(p["scale"]),
    shape    = as.numeric(p["shape"])
  ))
}

# Execute non-parametric bootstrap (2,000 iterations) on 28-day block maxima
set.seed(123)
boot_gev_res <- boot(data = final_maxima, statistic = boot_gev_lmom, R = 2000)

# Extract point estimates from the primary L-moments model fit
gev_pars <- distill(final_fit)

# Construct GEV parameter uncertainty summary table
gev_uncertainty_table <- data.frame(
  Parameter   = c("Location (μ)", "Scale (σ)", "Shape (ξ)"),
  Estimate    = round(as.numeric(gev_pars[c("location", "scale", "shape")]), 4),
  Std_Error   = round(apply(boot_gev_res$t, 2, sd, na.rm = TRUE), 4),
  Lower_95_CI = round(apply(boot_gev_res$t, 2, quantile, probs = 0.025, na.rm = TRUE), 4),
  Upper_95_CI = round(apply(boot_gev_res$t, 2, quantile, probs = 0.975, na.rm = TRUE), 4)
) %>%
  mutate(
    Formatted_Report = sprintf("%.4f (±%.4f)", Estimate, 1.96 * Std_Error),
    CI_95_Percent    = sprintf("[%.4f, %.4f]", Lower_95_CI, Upper_95_CI)
  )

# Display table in R console
cat("\n=== GEV L-MOMENTS PARAMETER UNCERTAINTY ===\n")
print(gev_uncertainty_table %>% select(Parameter, Estimate, Std_Error, CI_95_Percent), row.names = FALSE)

# Export LaTeX table
print(xtable(gev_uncertainty_table %>% select(Parameter, Estimate, Std_Error, CI_95_Percent),
             caption = "GEV parameter estimates ($\\\\mu, \\\\sigma, \\\\xi$) obtained via L-moments with 2,000 bootstrap standard errors and 95\\\\% confidence intervals.",
             label   = "tab:gev_lmom_uncertainty"),
      include.rownames = FALSE, booktabs = TRUE, comment = FALSE)







# ==========================================================
# 4. NON-STATIONARY GEV MODEL COMPARISON (WITH AICc)
# ==========================================================
# Prepare data and temporal covariates based on 28-day maxima
n_blocks <- length(final_maxima)
time_cov <- 1:n_blocks
time_cov_sq <- time_cov^2

# Fit candidate GEV models via Maximum Likelihood Estimation (MLE)
fit_stat <- fevd(final_maxima, type = "GEV")
fit_lin  <- fevd(final_maxima, type = "GEV", 
                 location.fun = ~time_cov, data = data.frame(time_cov))
fit_quad <- fevd(final_maxima, type = "GEV", 
                 location.fun = ~time_cov + time_cov_sq, 
                 data = data.frame(time_cov, time_cov_sq))

models <- list(Stationary = fit_stat, `Linear Trend` = fit_lin, `Quadratic Trend` = fit_quad)

# Extract Log-Likelihood, AIC, BIC, and Number of Parameters
loglik_vals <- sapply(models, function(x) -1 * as.numeric(x$results$value))
aic_vals    <- sapply(models, function(x) summary(x)$AIC)
bic_vals    <- sapply(models, function(x) summary(x)$BIC)
k_vals      <- c(Stationary = 3, `Linear Trend` = 4, `Quadratic Trend` = 5)

# Calculate Small-Sample Corrected Akaike Information Criterion (AICc)
aicc_vals   <- aic_vals + (2 * k_vals * (k_vals + 1)) / (n_blocks - k_vals - 1)

# Build unified model comparison dataframe
comp_results <- data.frame(
  Model = names(k_vals),
  `Log-Lik.` = round(loglik_vals, 3),
  AIC = round(aic_vals, 3),
  AICc = round(aicc_vals, 3),
  BIC = round(bic_vals, 3),
  `Par.` = k_vals,
  check.names = FALSE
)

# Display table in R console
cat("\n=== NON-STATIONARY GEV MODEL COMPARISON (AICc & BIC) ===\n")
print(comp_results, row.names = FALSE)

# Export LaTeX Table using xtable
print(xtable(comp_results, 
             digits = c(0, 0, 3, 3, 3, 3, 0), 
             caption = "Model comparison metrics for stationary and non-stationary GEV fits including log-likelihood, AIC, AICc, and BIC.",
             label = "tab:gev_nonstat_comparison"), 
      include.rownames = FALSE, 
      comment = FALSE,
      booktabs = TRUE)

# ==========================================================
# LIKELIHOOD-RATIO TESTS (LRT) FOR NESTED GEV MODELS
# ==========================================================

# 1. Extract Log-Likelihoods directly from fitted fevd objects
# Note: extRemes returns negative log-likelihoods in results$value
ll_stat <- -1 * as.numeric(fit_stat$results$value)
ll_lin  <- -1 * as.numeric(fit_lin$results$value)
ll_quad <- -1 * as.numeric(fit_quad$results$value)

# 2. Compute Likelihood Ratio Test Statistics (2 * Δll)
lr_stat_lin  <- 2 * (ll_lin - ll_stat)
lr_lin_quad  <- 2 * (ll_quad - ll_lin)
lr_stat_quad <- 2 * (ll_quad - ll_stat)

# 3. Compute p-values from Chi-Square Distribution
p_stat_lin  <- pchisq(lr_stat_lin,  df = 1, lower.tail = FALSE)
p_lin_quad  <- pchisq(lr_lin_quad,  df = 1, lower.tail = FALSE)
p_stat_quad <- pchisq(lr_stat_quad, df = 2, lower.tail = FALSE)

# 4. Construct Structured Output Table
lrt_results <- data.frame(
  Comparison = c(
    "Stationary vs Linear Trend",
    "Linear Trend vs Quadratic Trend",
    "Stationary vs Quadratic Trend"
  ),
  `df` = c(1, 1, 2),
  `LogLik_Null` = round(c(ll_stat, ll_lin, ll_stat), 3),
  `LogLik_Alt`  = round(c(ll_lin, ll_quad, ll_quad), 3),
  `LR_Stat`     = round(c(lr_stat_lin, lr_lin_quad, lr_stat_quad), 3),
  `p_value`     = c(p_stat_lin, p_lin_quad, p_stat_quad)
) %>%
  mutate(
    `Formatted_p` = ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value))
  )

# 5. Display Clean Console Output
cat("\n=== LIKELIHOOD-RATIO TESTS FOR NESTED GEV MODELS ===\n")
print(
  lrt_results %>% 
    select(Comparison, df, LogLik_Null, LogLik_Alt, LR_Stat, Formatted_p), 
  row.names = FALSE
)

# 6. Export Publication-Ready LaTeX Table
lrt_latex <- lrt_results %>%
  transmute(
    Comparison = Comparison,
    `$\\Delta \\text{df}$` = df,
    `$\\ell_{\\text{Null}}$` = sprintf("%.3f", LogLik_Null),
    `$\\ell_{\\text{Alt}}$`  = sprintf("%.3f", LogLik_Alt),
    `$\\Lambda$ (LR)`       = sprintf("%.3f", LR_Stat),
    `$p$-value`              = Formatted_p
  )

print(
  xtable(
    lrt_latex,
    caption = "Likelihood-ratio test statistics ($\\Lambda$) and associated $p$-values evaluating nested stationary and non-stationary GEV specifications.",
    label   = "tab:gev_lrt_tests",
    align   = c("l", "l", "c", "r", "r", "r", "r")
  ),
  include.rownames = FALSE,
  floating.environment = "table",
  table.placement = "H",
  size = "footnotesize",
  sanitize.text.function = identity,
  booktabs = TRUE,
  comment = FALSE
)


# ==========================================================
# GEV PROBABILITY (P-P) PLOT
# ==========================================================


# Extract fitted GEV parameters
gev_par <- distill(final_fit)

mu_hat    <- as.numeric(gev_par["location"])
sigma_hat <- as.numeric(gev_par["scale"])
xi_hat    <- as.numeric(gev_par["shape"])

# Sort observed block maxima
obs <- sort(final_maxima)

# Empirical probabilities
empirical <- ppoints(length(obs))

# Theoretical probabilities under the fitted GEV
theoretical <- pgev(
  obs,
  loc = mu_hat,
  scale = sigma_hat,
  shape = xi_hat
)

# Data frame
pp_data <- data.frame(
  Theoretical = theoretical,
  Empirical = empirical
)

# Probability Plot
gev_pp_plot <- ggplot(pp_data,
                      aes(x = Theoretical,
                          y = Empirical)) +
  geom_point(colour = "steelblue", size = 3) +
  geom_abline(intercept = 0,
              slope = 1,
              colour = "red",
              linetype = "dashed",
              linewidth = 0.8) +
  labs(
    title = "Probability Plot for the Fitted GEV Model",
    x = "Theoretical Probability",
    y = "Empirical Probability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(gev_pp_plot)






# ----------------------------------------------------------
#  PRIMARY STATIONARY GEV FIT (L-MOMENTS)
# ----------------------------------------------------------
blocks_per_year <- 365.25 / 28  # Conversion factor: B ≈ 13.0446

fit_lmom <- fevd(final_maxima, type = "GEV", method = "Lmoments")
pars_lmom <- distill(fit_lmom)

# ----------------------------------------------------------
# 2. BOOTSTRAP PARAMETER UNCERTAINTY (FOR L-MOMENTS)
# ----------------------------------------------------------
boot_gev_pars <- function(data, indices) {
  d <- data[indices]
  fit <- tryCatch(
    fevd(d, type = "GEV", method = "Lmoments"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(location = NA, scale = NA, shape = NA))
  p <- distill(fit)
  return(c(
    location = as.numeric(p["location"]), 
    scale    = as.numeric(p["scale"]), 
    shape    = as.numeric(p["shape"])
  ))
}

set.seed(123)
boot_pars_res <- boot(data = final_maxima, statistic = boot_gev_pars, R = 2000)

# Extract Standard Errors & 95% Percentile Confidence Intervals
se_loc   <- sd(boot_pars_res$t[, 1], na.rm = TRUE)
se_scale <- sd(boot_pars_res$t[, 2], na.rm = TRUE)
se_shape <- sd(boot_pars_res$t[, 3], na.rm = TRUE)

ci_loc   <- quantile(boot_pars_res$t[, 1], probs = c(0.025, 0.975), na.rm = TRUE)
ci_scale <- quantile(boot_pars_res$t[, 2], probs = c(0.025, 0.975), na.rm = TRUE)
ci_shape <- quantile(boot_pars_res$t[, 3], probs = c(0.025, 0.975), na.rm = TRUE)

# Parameter Summary Table with Uncertainty
lmom_param_table <- data.frame(
  Parameter   = c("Location (μ)", "Scale (σ)", "Shape (ξ)"),
  Estimate    = round(as.numeric(pars_lmom[c("location", "scale", "shape")]), 4),
  Std_Error   = round(c(se_loc, se_scale, se_shape), 4),
  Lower_95_CI = round(c(ci_loc[1], ci_scale[1], ci_shape[1]), 4),
  Upper_95_CI = round(c(ci_loc[2], ci_scale[2], ci_shape[2]), 4)
)

cat("\n=== GEV L-MOMENTS PARAMETER UNCERTAINTY ===\n")
print(lmom_param_table)

print(xtable(lmom_param_table, 
             caption = "GEV parameter estimates obtained via L-moments with bootstrap standard errors and 95\\% confidence intervals.",
             label = "tab:gev_lmom_uncertainty"), 
      include.rownames = FALSE, 
      booktabs = TRUE, 
      comment = FALSE)

# ----------------------------------------------------------
# 3. RECONCILING 28-DAY BLOCKS WITH ANNUAL RETURN LEVELS
# ----------------------------------------------------------
return_periods_years <- c(2, 5, 10, 20)
return_periods_blocks <- return_periods_years * blocks_per_year

# Bootstrap return levels using L-moments
boot_gev_rl <- function(data, indices, periods, B) {
  d <- data[indices]
  fit <- tryCatch(fevd(d, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, length(periods)))
  return(as.numeric(return.level(fit, return.period = periods * B)))
}

set.seed(123)
boot_rl_res <- boot(
  data = final_maxima, 
  statistic = boot_gev_rl, 
  R = 2000, 
  periods = return_periods_years, 
  B = blocks_per_year
)

reconciled_rl_df <- data.frame(
  Return_Period_Years = paste0(return_periods_years, "-year"),
  Block_Return_Period = round(return_periods_blocks, 2),
  Estimate            = round(boot_rl_res$t0, 2),
  Lower_95_CI         = round(apply(boot_rl_res$t, 2, quantile, probs = 0.025, na.rm = TRUE), 2),
  Upper_95_CI         = round(apply(boot_rl_res$t, 2, quantile, probs = 0.975, na.rm = TRUE), 2)
)

cat("\n=== RECONCILED ANNUAL RETURN LEVELS (L-MOMENTS) ===\n")
print(reconciled_rl_df)

print(xtable(reconciled_rl_df, 
             caption = "Stationary GEV return levels (cases per million) estimated via L-moments and reconciled to annual scale ($B = 13.0446$).",
             label = "tab:reconciled_lmom_rl"), 
      include.rownames = FALSE, 
      booktabs = TRUE, 
      comment = FALSE)
#===========================================================================================================

# ==========================================================

# ==========================================================

# ==========================================================
# MAIN GEV RETURN LEVEL PLOT & TABLE (TRUNCATED AT 20 YEARS)
# ==========================================================

# ==========================================================
# MAIN GEV RETURN LEVEL PLOT & TABLE (TRUNCATED AT 20 YEARS)
# ==========================================================

# 1. Annual conversion factor for 28-day blocks
blocks_per_year <- 365.25 / 28  # B ≈ 13.0446

# 2. Compute Gringorten empirical plotting positions
obs_sorted <- sort(final_maxima)
n_obs <- length(obs_sorted)
empirical_ranks <- 1:n_obs
empirical_periods_blocks <- (n_obs + 0.12) / (n_obs + 0.38 - empirical_ranks)
empirical_periods_years  <- empirical_periods_blocks / blocks_per_year

empirical_df <- data.frame(
  Years    = empirical_periods_years,
  Observed = obs_sorted
) %>%
  filter(Years <= 20)

# Determine minimum empirical return period for scale lower bound
min_year <- min(empirical_df$Years, na.rm = TRUE)

# 3. Primary discrete return periods & continuous evaluation grid
primary_years <- c(2, 5, 10, 20)
grid_periods_years <- exp(seq(log(min_year), log(20), length.out = 200))

# Combined bootstrap evaluation function
boot_gev_combined <- function(data, indices, grid_years, discrete_years, B) {
  d <- data[indices]
  fit <- tryCatch(fevd(d, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit)) {
    return(c(rep(NA, length(grid_years)), rep(NA, length(discrete_years))))
  }
  rl_grid <- as.numeric(return.level(fit, return.period = grid_years * B))
  rl_disc <- as.numeric(return.level(fit, return.period = discrete_years * B))
  return(c(rl_grid, rl_disc))
}

set.seed(123)
boot_combined_res <- boot(
  data = final_maxima, 
  statistic = boot_gev_combined, 
  R = 1000, 
  grid_years = grid_periods_years,
  discrete_years = primary_years,
  B = blocks_per_year
)

# Extract continuous ribbon data
n_grid <- length(grid_periods_years)
smooth_plot_df <- data.frame(
  Years        = grid_periods_years,
  Estimated_RL = boot_combined_res$t0[1:n_grid],
  Lower_CI     = pmax(0, apply(boot_combined_res$t[, 1:n_grid], 2, quantile, probs = 0.025, na.rm = TRUE)),
  Upper_CI     = apply(boot_combined_res$t[, 1:n_grid], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# Build the discrete summary table/points (replaces boot_summary_plot)
n_disc <- length(primary_years)
disc_indices <- (n_grid + 1):(n_grid + n_disc)

primary_summary_plot <- data.frame(
  Years        = primary_years,
  Estimated_RL = boot_combined_res$t0[disc_indices],
  Lower_CI     = pmax(0, apply(boot_combined_res$t[, disc_indices], 2, quantile, probs = 0.025, na.rm = TRUE)),
  Upper_CI     = apply(boot_combined_res$t[, disc_indices], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# ==========================================================
# PRINT NUMERIC TABLE (CONSOLE, MARKDOWN, LATEX)
# ==========================================================
primary_rl_table <- primary_summary_plot %>%
  transmute(
    `Return Period` = paste0(Years, "-year"),
    Estimate        = round(Estimated_RL, 2),
    `Lower 95% CI`  = round(Lower_CI, 2),
    `Upper 95% CI`  = round(Upper_CI, 2)
  )

cat("\n=== PRIMARY GEV RETURN LEVELS (2-20 YEARS) ===\n")
print(primary_rl_table, row.names = FALSE)

print(
  xtable(
    primary_rl_table,
    caption = "Primary stationary GEV return levels (cases per million) for 2-, 5-, 10-, and 20-year horizons.",
    label   = "tab:primary_gev_return_levels"
  ),
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE
)

# ==========================================================
# PUBLICATION-READY FIGURE
# ==========================================================
main_gev_rl_plot <- ggplot() +
  geom_ribbon(data = smooth_plot_df, 
              aes(x = Years, ymin = Lower_CI, ymax = Upper_CI), 
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = smooth_plot_df, 
            aes(x = Years, y = Lower_CI), 
            color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = smooth_plot_df, 
            aes(x = Years, y = Upper_CI), 
            color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = smooth_plot_df, 
            aes(x = Years, y = Estimated_RL), 
            color = "darkred", linewidth = 0.9) +
  geom_point(data = primary_summary_plot, 
             aes(x = Years, y = Estimated_RL), 
             color = "darkred", size = 2.5) +
  geom_point(data = empirical_df, 
             aes(x = Years, y = Observed), 
             shape = 21, fill = "white", color = "black", size = 2, stroke = 0.8) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20),
                labels = c("0.5", "1", "2", "5", "10", "20")) +
  coord_cartesian(xlim = c(min_year * 0.9, 22), ylim = c(0, NA)) +
  labs(
    x = "Return Period (Years)", 
    y = "Return Level (Cases per Million)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "plain")
  )

print(main_gev_rl_plot)

#==========================================================================================
#  Annual conversion factor for 28-day blocks
blocks_per_year <- 365.25 / 28  # B ≈ 13.0446

# 2. Compute Gringorten empirical plotting positions first
obs_sorted <- sort(final_maxima)
n_obs <- length(obs_sorted)
empirical_ranks <- 1:n_obs
empirical_periods_blocks <- (n_obs + 0.12) / (n_obs + 0.38 - empirical_ranks)
empirical_periods_years  <- empirical_periods_blocks / blocks_per_year

empirical_df <- data.frame(
  Years    = empirical_periods_years,
  Observed = obs_sorted
) %>%
  filter(Years <= 20)

# Determine the minimum empirical return period (between 0 and 1 year)
min_year <- min(empirical_df$Years, na.rm = TRUE)

# 3. Define fine evaluation grid starting from min_year up to 20 years
grid_periods_years <- exp(seq(log(min_year), log(20), length.out = 200))

# Bootstrap function for smooth confidence intervals
boot_smooth_fn <- function(data, indices, periods_grid, B) {
  d <- data[indices]
  fit <- tryCatch(fevd(d, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, length(periods_grid)))
  return(as.numeric(return.level(fit, return.period = periods_grid * B)))
}

set.seed(123)
boot_smooth_res <- boot(
  data = final_maxima, 
  statistic = boot_smooth_fn, 
  R = 1000, 
  periods_grid = grid_periods_years, 
  B = blocks_per_year
)

smooth_plot_df <- data.frame(
  Years        = grid_periods_years,
  Estimated_RL = boot_smooth_res$t0,
  Lower_CI     = pmax(0, apply(boot_smooth_res$t, 2, quantile, probs = 0.025, na.rm = TRUE)), # Constrained at 0
  Upper_CI     = apply(boot_smooth_res$t, 2, quantile, probs = 0.975, na.rm = TRUE)
)

# 4. Filter primary discrete points (2, 5, 10, 20 years)
primary_summary_plot <- boot_summary_plot %>%
  filter(Years <= 20)

# 5. Build Updated Figure (No Grid, Extended Steelblue Band below 1)
main_gev_rl_plot <- ggplot() +
  # Extended 95% Confidence Interval Band & Dashed Bounds
  geom_ribbon(data = smooth_plot_df, 
              aes(x = Years, ymin = Lower_CI, ymax = Upper_CI), 
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = smooth_plot_df, 
            aes(x = Years, y = Lower_CI), 
            color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = smooth_plot_df, 
            aes(x = Years, y = Upper_CI), 
            color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  
  # Smooth Fitted GEV Return Level Model Line
  geom_line(data = smooth_plot_df, 
            aes(x = Years, y = Estimated_RL), 
            color = "darkred", linewidth = 0.9) +
  
  # Target Primary Return Period Points
  geom_point(data = primary_summary_plot, 
             aes(x = Years, y = Estimated_RL), 
             color = "darkred", size = 2.5) +
  
  # All Empirical Points Overlaid (Including those between 0 and 1)
  geom_point(data = empirical_df, 
             aes(x = Years, y = Observed), 
             shape = 21, fill = "white", color = "black", size = 2, stroke = 0.8) +
  
  # Log Scale Display
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20),
                labels = c("0.5", "1", "2", "5", "10", "20")) +
  coord_cartesian(xlim = c(min_year * 0.9, 22), ylim = c(0, NA)) +
  
  labs(
    x = "Return Period (Years)", 
    y = "Return Level (Cases per Million)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(), # Grid lines removed
    axis.title = element_text(face = "plain")
  )

print(main_gev_rl_plot)



# Ensure required vectors/dataframes exist
# Expected input: `zambia_data` with columns `Day` and `cases`

# ------------------------------------------------------------------------------
# . DATA PREPROCESSING & TIME HORIZON SETUP
# ------------------------------------------------------------------------------
# Extract positive incidence values for upper-tail analysis
active_cases <- zambia_data %>%
  filter(!is.na(cases), cases > 0) %>%
  pull(cases)

n_positive_obs <- length(active_cases)

# Calculate exact temporal span of the study period in years
study_days <- as.numeric(
  max(zambia_data$Day, na.rm = TRUE) - min(zambia_data$Day, na.rm = TRUE)
) + 1

study_years <- study_days / 365.25

# Define candidate quantiles (85th to 95th percentiles)
percentiles <- seq(0.85, 0.95, by = 0.01)
u_candidates <- quantile(active_cases, probs = percentiles, na.rm = TRUE)

# ------------------------------------------------------------------------------
# 3. THRESHOLD DIAGNOSTIC LOOP
# ------------------------------------------------------------------------------
evaluate_threshold <- function(u, percentile_label, data, n_total, years, run_length = 7) {
  # Raw exceedances
  exceedances <- data[data > u]
  raw_N       <- length(exceedances)
  
  if (raw_N == 0) {
    return(data.frame(
      Percentile           = percentile_label,
      Threshold_u          = round(u, 2),
      Raw_N                = 0,
      Clusters_Nc          = 0,
      Theta                = NA_real_,
      Raw_Rate_Daily       = 0,
      Cluster_Rate_Annual  = 0,
      Scale_sigma          = NA_real_,
      Shape_xi             = NA_real_
    ))
  }
  
  # Temporal Declustering (Runs Method, r = 7 days)
  declustered_series <- decluster(
    x         = data,
    threshold = u,
    method    = "runs",
    r         = run_length
  )
  
  # Count independent cluster peaks
  cluster_N <- sum(declustered_series > u, na.rm = TRUE)
  
  # Extremal Index (Theta = Reciprocal of Mean Cluster Size)
  theta <- cluster_N / raw_N
  
  # Exceedance Rates
  raw_rate_daily  <- raw_N / n_total
  lambda_c_annual <- cluster_N / years
  
  # Fit Generalized Pareto Distribution (GPD) via L-moments
  gpd_fit <- tryCatch({
    fevd(
      x         = declustered_series,
      threshold = u,
      type      = "GP",
      method    = "Lmoments"
    )
  }, error = function(e) NULL)
  
  # Parameter Extraction
  if (is.null(gpd_fit)) {
    sigma_hat <- NA_real_
    xi_hat    <- NA_real_
  } else {
    pars      <- distill(gpd_fit)
    sigma_hat <- as.numeric(pars["scale"])
    xi_hat    <- as.numeric(pars["shape"])
  }
  
  # Construct Output Row
  data.frame(
    Percentile           = percentile_label,
    Threshold_u          = round(u, 2),
    Raw_N                = raw_N,
    Clusters_Nc          = cluster_N,
    Theta                = round(theta, 3),
    Raw_Rate_Daily       = round(raw_rate_daily, 4),
    Cluster_Rate_Annual  = round(lambda_c_annual, 4),
    Scale_sigma          = round(sigma_hat, 4),
    Shape_xi             = round(xi_hat, 4)
  )
}

# Run diagnostics across all candidate thresholds
threshold_diagnostics_df <- map2_df(
  .x = u_candidates,
  .y = paste0(round(percentiles * 100), "%"),
  .f = ~ evaluate_threshold(
    u                = .x,
    percentile_label = .y,
    data             = active_cases,
    n_total          = n_positive_obs,
    years            = study_years,
    run_length       = 7
  )
)

# ------------------------------------------------------------------------------
# 4. CONSOLE & PUBLICATION-READY OUTPUTS
# ------------------------------------------------------------------------------
# Clean Console Output
cat("\n=== POT THRESHOLD SELECTION DIAGNOSTICS ===\n")
print(threshold_diagnostics_df, row.names = FALSE)

# Publication-Ready Markdown Table
kable(
  threshold_diagnostics_df,
  caption   = "Threshold Selection Diagnostics and GPD Parameter Estimates (L-Moments)",
  col.names = c("Percentile", "Threshold (u)", "Raw N", "Clusters (Nc)", 
                "Theta (θ)", "Daily Raw Rate", "Annual Cluster Rate (λc)", 
                "Scale (σ)", "Shape (ξ)"),
  align     = c("l", "r", "r", "r", "r", "r", "r", "r", "r")
)

# Publication-Ready LaTeX Output
latex_table <- xtable(
  threshold_diagnostics_df,
  caption = "Threshold diagnostics and Generalised Pareto Distribution parameter estimates across candidate percentiles.",
  label   = "tab:pot_threshold_diagnostics",
  digits  = c(0, 0, 2, 0, 0, 3, 4, 4, 4, 4)
)

# Print LaTeX Table
print(
  latex_table,
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE
)



# =========================================================================================
#  SELECTED POT THRESHOLD & CLUSTER ANALYSIS
# ==============================================================================

# Parameters
u_selected <- 61.26
run_length <- 7

# ------------------------------------------------------------------------------
# 1. Identify and sequence threshold exceedances
# ------------------------------------------------------------------------------
exceedance_data <- zambia_data %>%
  filter(cases > u_selected) %>%
  mutate(Day = as.Date(Day)) %>%
  arrange(Day) %>%
  mutate(
    gap_days    = as.numeric(Day - lag(Day)),
    new_cluster = if_else(is.na(gap_days) | gap_days > run_length, 1, 0),
    cluster_id  = cumsum(new_cluster)
  )

# Guard check: proceed only if exceedances exist
if (nrow(exceedance_data) == 0) {
  
  cat("No exceedances found above threshold:", u_selected, "\n")
  
} else {
  
  # ----------------------------------------------------------------------------
  # 2. Construct cluster summary
  # ----------------------------------------------------------------------------
  cluster_data <- exceedance_data %>%
    group_by(cluster_id) %>%
    summarise(
      Cluster_Start   = min(Day),
      Cluster_End     = max(Day),
      Duration_Days   = as.numeric(Cluster_End - Cluster_Start) + 1,
      Cluster_Size    = n(),
      Cluster_Maximum = max(cases),
      .groups         = "drop"
    )
  
  # ----------------------------------------------------------------------------
  # 3. Calculate summary statistics
  # ----------------------------------------------------------------------------
  raw_N_selected       <- nrow(exceedance_data)
  cluster_N_selected   <- nrow(cluster_data)
  theta_selected       <- cluster_N_selected / raw_N_selected
  mean_cluster_duration <- mean(cluster_data$Duration_Days)
  mean_cluster_size     <- mean(cluster_data$Cluster_Size)
  
  # ----------------------------------------------------------------------------
  # 4. Console output
  # ----------------------------------------------------------------------------
  cat(
    "\n==========================================================\n",
    "SELECTED POT CLUSTER ANALYSIS\n",
    "==========================================================\n",
    "Selected threshold:   ", u_selected, "\n",
    "Runs parameter:       ", run_length, " days\n",
    "Raw exceedances:      ", raw_N_selected, "\n",
    "Independent clusters: ", cluster_N_selected, "\n",
    "Extremal index:       ", round(theta_selected, 4), "\n",
    "Mean cluster duration:", round(mean_cluster_duration, 2), " days\n",
    "Mean cluster size:    ", round(mean_cluster_size, 2), " exceedances\n",
    "==========================================================\n\n",
    "Cluster details:\n",
    sep = ""
  )
  
  print(cluster_data, row.names = FALSE)
}



# ===================================================================================================
#  RUN-LENGTH SENSITIVITY ANALYSIS
# ==============================================================================

# Define candidate run lengths (days)
run_lengths <- c(1, 2, 3, 4, 5, 7, 10, 14, 21)

# Evaluate declustering across run lengths
run_length_sensitivity <- map_df(run_lengths, function(r) {
  
  # Perform declustering using the runs method
  dec <- decluster(
    x         = active_cases,
    threshold = u_selected,
    method    = "runs",
    r         = r
  )
  
  # Calculate cluster metrics
  cluster_N <- sum(dec > u_selected, na.rm = TRUE)
  theta     <- cluster_N / raw_N_selected
  lambda_c  <- cluster_N / study_years
  
  # Assemble output row
  tibble(
    Run_Length_Days      = r,
    Raw_Exceedances      = raw_N_selected,
    Independent_Clusters = cluster_N,
    Theta                = round(theta, 4),
    Annual_Cluster_Rate  = round(lambda_c, 4)
  )
})

# ------------------------------------------------------------------------------
# Output results
# ------------------------------------------------------------------------------
cat(
  "\n==========================================================\n",
  "RUN-LENGTH SENSITIVITY ANALYSIS\n",
  "==========================================================\n",
  sep = ""
)

print(run_length_sensitivity, row.names = FALSE)


#===========================================================================================================
# ==============================================================================
# Peak-Over-Threshold (POT) & Extremal Index Analysis
# ==============================================================================

# Required Libraries


# ------------------------------------------------------------------------------
# 1. Configuration & Threshold Settings
# ------------------------------------------------------------------------------

u_selected <- 61.26
run_length <- 7

# ------------------------------------------------------------------------------
# 2. Identify Threshold Exceedances & Decluster into Clusters
# ------------------------------------------------------------------------------

exceedance_data <- zambia_data %>%
  filter(cases > u_selected) %>%
  arrange(Day) %>%
  mutate(
    gap_days    = as.numeric(Day - lag(Day)),
    new_cluster = ifelse(is.na(gap_days) | gap_days > run_length, 1, 0),
    cluster_id  = cumsum(new_cluster)
  )

# Summarize Cluster Characteristics
cluster_data <- exceedance_data %>%
  group_by(cluster_id) %>%
  summarise(
    Cluster_Start   = min(Day),
    Cluster_End     = max(Day),
    Duration_Days   = as.numeric(Cluster_End - Cluster_Start) + 1,
    Cluster_Size    = n(),
    Cluster_Maximum = max(cases),
    .groups         = "drop"
  )

# ------------------------------------------------------------------------------
# 3. Compute Summary Cluster Statistics
# ------------------------------------------------------------------------------

raw_N_selected        <- nrow(exceedance_data)
cluster_N_selected    <- nrow(cluster_data)
theta_selected        <- cluster_N_selected / raw_N_selected
mean_cluster_duration <- mean(cluster_data$Duration_Days)
mean_cluster_size     <- mean(cluster_data$Cluster_Size)

cat("\n==========================================================\n")
cat("SELECTED POT CLUSTER ANALYSIS\n")
cat("==========================================================\n")
cat("Selected threshold:   ", u_selected, "\n")
cat("Runs parameter:       ", run_length, "days\n")
cat("Raw exceedances:      ", raw_N_selected, "\n")
cat("Independent clusters: ", cluster_N_selected, "\n")
cat("Extremal index (θ):   ", round(theta_selected, 4), "\n")
cat("Mean cluster duration:", round(mean_cluster_duration, 2), "days\n")
cat("Mean cluster size:    ", round(mean_cluster_size, 2), "exceedances\n\n")

cat("Cluster Details:\n")
print(cluster_data, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Generalized Pareto Distribution (GPD) Fit
# ------------------------------------------------------------------------------

cluster_peaks <- cluster_data$Cluster_Maximum

final_gpd_fit <- fevd(
  x         = cluster_peaks,
  threshold = u_selected,
  type      = "GP",
  method    = "Lmoments"
)

gpd_parameters <- distill(final_gpd_fit)
sigma_hat      <- as.numeric(gpd_parameters["scale"])
xi_hat         <- as.numeric(gpd_parameters["shape"])

cat("\n==========================================================\n")
cat("SELECTED GPD PARAMETERS\n")
cat("==========================================================\n")
cat("Scale (σ):", round(sigma_hat, 4), "\n")
cat("Shape (ξ):", round(xi_hat, 4), "\n")

# ------------------------------------------------------------------------------
# 5. Run-Length Sensitivity Analysis
# ------------------------------------------------------------------------------

run_lengths <- c(1, 2, 3, 4, 5, 7, 10, 14, 21)

run_length_sensitivity <- map_df(run_lengths, function(r) {
  dec <- decluster(
    x         = active_cases,
    threshold = u_selected,
    method    = "runs",
    r         = r
  )
  
  cluster_N <- sum(dec > u_selected, na.rm = TRUE)
  theta     <- cluster_N / raw_N_selected
  lambda_c  <- cluster_N / study_years
  
  data.frame(
    Run_Length_Days      = r,
    Raw_Exceedances      = raw_N_selected,
    Independent_Clusters = cluster_N,
    Theta                = round(theta, 4),
    Annual_Cluster_Rate  = round(lambda_c, 4)
  )
})

cat("\n==========================================================\n")
cat("RUN-LENGTH SENSITIVITY\n")
cat("==========================================================\n")
print(run_length_sensitivity, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Extremal Index Uncertainty (Nonparametric Bootstrap)
# ------------------------------------------------------------------------------

set.seed(123)

cluster_sizes <- cluster_data$Cluster_Size
B_theta       <- 2000

theta_boot <- replicate(B_theta, {
  sampled_sizes <- sample(cluster_sizes, size = length(cluster_sizes), replace = TRUE)
  length(sampled_sizes) / sum(sampled_sizes)
})

theta_CI <- quantile(theta_boot, probs = c(0.025, 0.975), na.rm = TRUE)

cat("\n==========================================================\n")
cat("EXTREMAL INDEX UNCERTAINTY\n")
cat("==========================================================\n")
cat("Theta Estimate:", round(theta_selected, 4), "\n")
cat("Lower 95% CI:  ", round(theta_CI[1], 4), "\n")
cat("Upper 95% CI:  ", round(theta_CI[2], 4), "\n")

# ------------------------------------------------------------------------------
# 7. GPD Parameter Uncertainty (Nonparametric Bootstrap)
# ------------------------------------------------------------------------------

set.seed(123)

B_gpd <- 2000

gpd_boot <- map_dfr(1:B_gpd, function(i) {
  sampled_peaks <- sample(cluster_peaks, size = length(cluster_peaks), replace = TRUE)
  
  fit_boot <- tryCatch(
    fevd(
      x         = sampled_peaks,
      threshold = u_selected,
      type      = "GP",
      method    = "Lmoments"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_boot)) {
    return(data.frame(Sigma = NA_real_, Xi = NA_real_))
  }
  
  pars_boot <- distill(fit_boot)
  
  data.frame(
    Sigma = as.numeric(pars_boot["scale"]),
    Xi    = as.numeric(pars_boot["shape"])
  )
})

# Filter finite values for percentile interval computation
finite_sigma <- gpd_boot$Sigma[is.finite(gpd_boot$Sigma)]
finite_xi    <- gpd_boot$Xi[is.finite(gpd_boot$Xi)]

sigma_CI <- quantile(finite_sigma, probs = c(0.025, 0.975), na.rm = TRUE)
xi_CI    <- quantile(finite_xi, probs = c(0.025, 0.975), na.rm = TRUE)

# Combine parameter estimates into a summary matrix
parameter_uncertainty <- data.frame(
  Parameter   = c("Theta (θ)", "Scale (σ)", "Shape (ξ)"),
  Estimate    = round(c(theta_selected, sigma_hat, xi_hat), 4),
  Lower_95_CI = round(c(theta_CI[1], sigma_CI[1], xi_CI[1]), 4),
  Upper_95_CI = round(c(theta_CI[2], sigma_CI[2], xi_CI[2]), 4)
)

cat("\n==========================================================\n")
cat("POT PARAMETER UNCERTAINTY\n")
cat("==========================================================\n")
print(parameter_uncertainty, row.names = FALSE)
cat("\nFinite bootstrap sigma estimates:", length(finite_sigma), "of", B_gpd, "\n")
cat("Finite bootstrap xi estimates:   ", length(finite_xi), "of", B_gpd, "\n")

# ------------------------------------------------------------------------------
# 8. Overall POT Rate Summary
# ------------------------------------------------------------------------------

lambda_cluster_annual <- cluster_N_selected / study_years
raw_rate_daily        <- raw_N_selected / n_positive_obs

pot_summary <- data.frame(
  Threshold                  = u_selected,
  Run_Length_Days            = run_length,
  Raw_Exceedances            = raw_N_selected,
  Independent_Clusters      = cluster_N_selected,
  Theta                      = round(theta_selected, 4),
  Raw_Rate_Daily             = round(raw_rate_daily, 4),
  Cluster_Rate_Annual        = round(lambda_cluster_annual, 4),
  Sigma                      = round(sigma_hat, 4),
  Xi                         = round(xi_hat, 4),
  Mean_Cluster_Duration_Days = round(mean_cluster_duration, 2),
  Mean_Cluster_Size          = round(mean_cluster_size, 2)
)

cat("\n==========================================================\n")
cat("FINAL POT SUMMARY\n")
cat("==========================================================\n")
print(pot_summary, row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. Publication-Ready LaTeX Tables
# ------------------------------------------------------------------------------

cat("\n==========================================================\n")
cat("LATEX: CLUSTER TABLE\n")
cat("==========================================================\n")
print(
  xtable(
    cluster_data,
    caption = "Identified exceedance clusters at the selected threshold using a 7-day runs parameter.",
    label   = "tab:cluster_details"
  ),
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE,
  digits           = 2
)

cat("\n==========================================================\n")
cat("LATEX: RUN-LENGTH SENSITIVITY TABLE\n")
cat("==========================================================\n")
print(
  xtable(
    run_length_sensitivity,
    caption = "Sensitivity of runs declustering to alternative run lengths.",
    label   = "tab:run_length_sensitivity"
  ),
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE,
  digits           = 4
)

cat("\n==========================================================\n")
cat("LATEX: PARAMETER UNCERTAINTY TABLE\n")
cat("==========================================================\n")
print(
  xtable(
    parameter_uncertainty,
    caption = "Bootstrap uncertainty estimates for the extremal index and GPD parameters at the selected threshold.",
    label   = "tab:pot_uncertainty"
  ),
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE,
  digits           = 4
)

cat("\n==========================================================\n")
cat("LATEX: PRIMARY RETURN LEVEL TABLE\n")
cat("==========================================================\n")
print(
  xtable(
    gev_return_table,
    caption = "Primary stationary GEV return levels for 2-, 5-, 10- and 20-year return periods.",
    label   = "tab:gev_return_levels"
  ),
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE,
  digits           = 2
)

cat("\n==========================================================\n")
cat("LATEX: LONG-HORIZON SENSITIVITY TABLE\n")
cat("==========================================================\n")
print(
  xtable(
    long_return_table,
    caption = "Long-horizon stationary GEV return levels presented as sensitivity analysis.",
    label   = "tab:long_horizon_return_levels"
  ),
  include.rownames = FALSE,
  booktabs         = TRUE,
  comment          = FALSE,
  digits           = 2
)

# ------------------------------------------------------------------------------
# 10. Final Consistency Verification Check
# ------------------------------------------------------------------------------

cat("\n==========================================================\n")
cat("FINAL ANALYSIS CHECK\n")
cat("==========================================================\n")
cat("Primary block size:     28 days\n")
cat("Blocks per year:        ", round(blocks_per_year, 4), "\n")
cat("Selected POT threshold: ", u_selected, "\n")
cat("Runs parameter:         ", run_length, "days\n")
cat("Raw exceedances:        ", raw_N_selected, "\n")
cat("Independent clusters:   ", cluster_N_selected, "\n")
cat("Theta (θ):              ", round(theta_selected, 4), "\n")
cat("Mean cluster duration:  ", round(mean_cluster_duration, 2), "days\n")
cat("Mean cluster size:      ", round(mean_cluster_size, 2), "\n")
cat("GPD Scale (σ):          ", round(sigma_hat, 4), "\n")
cat("GPD Shape (ξ):          ", round(xi_hat, 4), "\n")
cat("\nAnalysis completed successfully.\n")


# ----------------------------------------------------------
# We extract dates and cases directly from the filtered dataframe in Section 1
dates <- zambia_data$Day
active_cases <- zambia_data$cases

u_selected <- 61.26

# ----------------------------------------------------------
#  UPDATED TIME SERIES EXCEEDANCES PLOT
# ----------------------------------------------------------
# Prepare the plotting dataframe with the synchronized dates
plot_ts_df <- data.frame(Date = dates, Cases = active_cases) %>%
  mutate(Status = ifelse(Cases > u_selected, "Exceedance", "Baseline"))

ts_gg <- ggplot(plot_ts_df, aes(x = Date, y = Cases)) +
  # Draw the base line for all cases
  geom_line(color = "gray80", linewidth = 0.5) +
  # Highlight the exceedances in dark red
  geom_point(data = filter(plot_ts_df, Status == "Exceedance"), 
             aes(color = Status), size = 1.5) +
  # The threshold line
  geom_hline(yintercept = u_selected, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("Exceedance" = "darkred")) +
  # FORCING THE X-AXIS TO SHOW 2020 THROUGH 2023
  scale_x_date(
    date_breaks = "1 year", 
    date_labels = "%Y",
    limits = c(as.Date("2020-01-01"), as.Date("2023-12-31"))
  ) +
  labs(
    title = "COVID-19 Exceedances in Zambia",
    x = "Timeline (2020 - 2023)", 
    y = "Active Cases (per Million)"
  ) +
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Print the fixed plot
print(ts_gg)





#  Time Series Exceedances plot
plot_ts_df <- data.frame(Date = dates, Cases = active_cases) %>%
  mutate(Status = ifelse(Cases > u_selected, "Exceedance", "Baseline"))

ts_gg <- ggplot(plot_ts_df, aes(x = Date, y = Cases)) +
  geom_line(color = "gray80", linewidth = 0.5) +
  geom_point(data = filter(plot_ts_df, Status == "Exceedance"), 
             aes(color = Status), size = 1.5) +
  geom_hline(yintercept = u_selected, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("Exceedance" = "darkred")) +
  labs(x = "Timeline", y = "Active Cases (per Million)") +
  theme_bw() + 
  theme(legend.position = "none", panel.grid = element_blank())

#  Mean Residual Life Plot
mrl_gg <- ggplot(threshold_table_full, aes(x = Threshold_u, y = Mean_Excess)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(aes(size = Clusters_Nc), color = "darkred") +
  geom_vline(xintercept = u_selected, linetype = "dashed") +
  labs(x = "Threshold (u)", y = "Mean Excess") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

# FIG 4.2: Shape Parameter Stability
shape_stab_gg <- ggplot(threshold_table_full, aes(x = Threshold_u, y = Shape_xi)) +
  geom_ribbon(aes(ymin = Shape_xi - 0.1, ymax = Shape_xi + 0.1), fill = "darkgreen", alpha = 0.1) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2) +
  geom_vline(xintercept = u_selected, linetype = "dashed", color = "red") +
  labs(x = "Threshold (u)", y = "Shape Estimate (\u03BE)") +
  theme_bw() + theme(panel.grid = element_blank())

#  Scale Parameter Stability
scale_stab_gg <- ggplot(threshold_table_full, aes(x = Threshold_u, y = Scale_sigma)) +
  geom_ribbon(aes(ymin = Scale_sigma * 0.9, ymax = Scale_sigma * 1.1), fill = "darkblue", alpha = 0.1) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 2) +
  geom_vline(xintercept = u_selected, linetype = "dashed", color = "red") +
  labs(x = "Threshold (u)", y = "Scale Estimate (\u03C3)") +
  theme_bw() + theme(panel.grid = element_blank())

# 100-Year Projection (Ref: gev_return.pdf)
plot_years <- seq(1, 100, length.out = 100)
rl_curve <- return.level(final_gpd_fit, return.period = plot_years * 365.25, do.ci = TRUE)
df_100yr <- data.frame(Year = plot_years, Level = rl_curve[,2], Lower = rl_curve[,1], Upper = rl_curve[,3])

rl_100yr_gg <- ggplot(df_100yr, aes(x = Year, y = Level)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkblue", alpha = 0.15) +
  geom_line(color = "darkblue", linewidth = 1) +
  labs(x = "Return Period (Years)", y = "Return Level (Cases per Million)") +
  theme_bw() + theme(panel.grid = element_blank())

#  Q-Q Plot (Ref: Q-Q PLOT.pdf)
all_data <- as.numeric(datagrabber(final_gpd_fit)[, 1])
exceedances <- sort(all_data[all_data > u_selected] - u_selected)
n_exc <- length(exceedances)
p_points <- (1:n_exc) / (n_exc + 1)
fit_p <- distill(final_gpd_fit)
theo_q <- (as.numeric(fit_p["scale"]) / as.numeric(fit_p["shape"])) * 
  ((1 - p_points)^(-as.numeric(fit_p["shape"])) - 1)

qq_gg <- ggplot(data.frame(Emp = exceedances, Theo = sort(theo_q)), aes(x = Theo, y = Emp)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(color = "steelblue", size = 3) +
  labs(x = "Theoretical Quantiles (GPD)", y = "Empirical Quantiles") +
  theme_bw() + theme(panel.grid = element_blank())

# ----------------------------------------------------------
#  PRINT ALL FIGURES
# ----------------------------------------------------------
print(ts_gg)
print(mrl_gg)
print(shape_stab_gg)
print(scale_stab_gg)
print(qq_gg)

# ==========================================================
# GPD PROBABILITY (P-P) PLOT
# ==========================================================

library(extRemes)
library(ggplot2)

# Extract fitted parameters
gpd_par <- distill(final_gpd_fit)

sigma_hat <- as.numeric(gpd_par["scale"])
xi_hat    <- as.numeric(gpd_par["shape"])

# Exceedances above the threshold
exceedances <- declustered_data[declustered_data > u_selected] - u_selected

# Empirical probabilities
n <- length(exceedances)
empirical <- ppoints(n)

# Theoretical GPD probabilities
theoretical <- pgpd(sort(exceedances),
                    loc = 0,
                    scale = sigma_hat,
                    shape = xi_hat)

# Data frame
pp_data <- data.frame(
  Empirical = empirical,
  Theoretical = theoretical
)

# P-P Plot
pp_plot <- ggplot(pp_data,
                  aes(x = Theoretical,
                      y = Empirical)) +
  geom_point(size = 3, colour = "steelblue") +
  geom_abline(intercept = 0,
              slope = 1,
              colour = "red",
              linetype = "dashed",
              linewidth = 0.8) +
  labs(
    title = "Probability Plot for the Fitted GPD Model",
    x = "Theoretical Probability",
    y = "Empirical Probability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(pp_plot)


# ------------------------------------------------------------------------------
# 8. GPD RETURN LEVEL BOOTSTRAP (95% CI)
# ------------------------------------------------------------------------------

# Return periods for the primary POT return-level analysis
return_periods <- c(2, 5, 10, 20)

# Annual cluster rate from the selected POT analysis
lambda_c <- lambda_cluster_annual

# GPD return-level calculation
calc_gpd_rl <- function(T, u, sigma, xi, rate) {
  
  if (abs(xi) < 1e-8) {
    return(u + sigma * log(rate * T))
  } else {
    return(
      u + (sigma / xi) * ((rate * T)^xi - 1)
    )
  }
}

# Point estimates of GPD return levels
estimated_RL <- sapply(
  return_periods,
  calc_gpd_rl,
  u = u_selected,
  sigma = sigma_hat,
  xi = xi_hat,
  rate = lambda_c
)

# ------------------------------------------------------------------------------
# Robust bootstrap function
# ------------------------------------------------------------------------------

gpd_boot_stat <- function(data, indices) {
  
  sample_peaks <- data[indices]
  
  # Ensure sufficient unique values for parameter estimation
  if (length(unique(sample_peaks)) < 3) {
    return(rep(NA_real_, length(return_periods)))
  }
  
  # Fit GPD to resampled cluster peaks
  fit_b <- tryCatch({
    fevd(
      x         = sample_peaks,
      threshold = u_selected,
      type      = "GP",
      method    = "Lmoments"
    )
  }, error = function(e) NULL)
  
  if (is.null(fit_b)) {
    return(rep(NA_real_, length(return_periods)))
  }
  
  # Extract bootstrap parameters
  pars_b <- tryCatch(
    distill(fit_b),
    error = function(e) NULL
  )
  
  if (is.null(pars_b)) {
    return(rep(NA_real_, length(return_periods)))
  }
  
  sigma_b <- as.numeric(pars_b["scale"])
  xi_b    <- as.numeric(pars_b["shape"])
  
  if (is.na(sigma_b) || is.na(xi_b) || sigma_b <= 0) {
    return(rep(NA_real_, length(return_periods)))
  }
  
  # Calculate return levels using bootstrap GPD parameters
  sapply(
    return_periods,
    calc_gpd_rl,
    u     = u_selected,
    sigma = sigma_b,
    xi    = xi_b,
    rate  = lambda_c
  )
}

# ------------------------------------------------------------------------------
# Execute bootstrap simulation
# ------------------------------------------------------------------------------

set.seed(123)

boot_results <- boot(
  data      = cluster_peaks,
  statistic = gpd_boot_stat,
  R         = 5000
)

# ------------------------------------------------------------------------------
# Extract valid bootstrap iterations
# ------------------------------------------------------------------------------

boot_matrix <- boot_results$t

valid_rows <- complete.cases(boot_matrix)

valid_boot <- boot_matrix[
  valid_rows,
  ,
  drop = FALSE
]

cat(
  "Successful Bootstrap Replicates:",
  sum(valid_rows),
  "out of",
  nrow(boot_matrix),
  "\n\n"
)

# ------------------------------------------------------------------------------
# Calculate 95% percentile confidence intervals
# ------------------------------------------------------------------------------

lower_CI <- apply(
  valid_boot,
  2,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

upper_CI <- apply(
  valid_boot,
  2,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

# ------------------------------------------------------------------------------
# Publication-ready return-level table
# ------------------------------------------------------------------------------

return_level_table <- data.frame(
  Return_Period = paste0(return_periods, "-year"),
  Lower_95_CI   = round(lower_CI, 2),
  Estimate      = round(estimated_RL, 2),
  Upper_95_CI   = round(upper_CI, 2)
)

# Display Markdown table
kable(
  return_level_table,
  caption = "POT-GPD Annual Return Level Estimates with 95% Bootstrap CIs",
  col.names = c(
    "Return Period",
    "Lower 95% CI",
    "Estimate",
    "Upper 95% CI"
  ),
  align = c("l", "r", "r", "r")
)

# ------------------------------------------------------------------------------
# LaTeX table
# ------------------------------------------------------------------------------

latex_rl_table <- xtable(
  return_level_table,
  caption = "POT-GPD Return level estimates (cases per million) with 95\\% confidence intervals.",
  label   = "tab:gpd_return_levels",
  digits  = c(0, 0, 2, 2, 2)
)

print(
  latex_rl_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  comment = FALSE
)



# ==========================================================
# 1. SETUP & DATA LOADING
# ==========================================================
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

# Data Preparation
# Filter and construct incidence vector for Zambia
zambia_data <- Cov_Za %>%
  mutate(Day = as.Date(Day),
         cases = as.numeric(Daily.new.confirmed.cases.of.COVID.19.per.million.people)) %>%
  filter(Entity == "Zambia",
         Day >= "2020-01-04" & Day <= "2023-12-31",
         !is.na(cases)) %>%
  arrange(Day)

cases_vector <- zambia_data$cases

# ==========================================================
# 2. ANALYSIS LOOP (BLOCK MAXIMA SELECTION)
# ==========================================================
results <- data.frame()
block_sizes <- c(7, 14, 21, 28)

for (b in block_sizes) {
  n_full_blocks <- floor(length(cases_vector) / b)
  if (n_full_blocks < 5) next
  
  trimmed_data <- cases_vector[1:(n_full_blocks * b)]
  data_matrix <- matrix(trimmed_data, nrow = b, byrow = FALSE)
  block_maxima_temp <- apply(data_matrix, 2, max, na.rm = TRUE)
  block_maxima_temp <- block_maxima_temp[block_maxima_temp > 0]
  
  tryCatch({
    # Fit GEV model using L-moments
    gev_fit <- fevd(block_maxima_temp, type = "GEV", method = "Lmoments")
    params <- distill(gev_fit)
    
    L <- as.numeric(params["location"])
    S <- as.numeric(params["scale"])
    G <- as.numeric(params["shape"])
    
    # Validation Statistics
    p_seq <- ppoints(length(block_maxima_temp))
    q_theory <- qgev(p_seq, loc = L, scale = S, shape = G)
    
    # RMSE
    rmse_val <- sqrt(mean((sort(block_maxima_temp) - q_theory)^2))
    
    # KS Test with jitter for ties
    jittered_maxima <- jitter(block_maxima_temp, amount = 1e-5)
    ks_pval <- ks.test(jittered_maxima, "pgev", loc = L, scale = S, shape = G)$p.value
    
    # Anderson-Darling Test
    ad_pval <- goftest::ad.test(block_maxima_temp, "pgev", loc = L, scale = S, shape = G)$p.value
    
    # Runs Test for Independence
    diffs <- diff(block_maxima_temp)
    signs <- sign(diffs)[sign(diffs) != 0]
    runs_pval <- if(length(signs) > 5) {
      tseries::runs.test(as.factor(signs > 0))$p.value
    } else {
      NA
    }
    
    # Wilcoxon Signed-Rank Test
    theo_median <- qgev(0.5, loc = L, scale = S, shape = G)
    wilcox_pval <- wilcox.test(block_maxima_temp, mu = theo_median)$p.value
    
    results <- rbind(results, data.frame(
      BlockSize = b, 
      nBlocks   = length(block_maxima_temp), 
      RMSE      = rmse_val,
      Runs_p    = runs_pval,
      Wilcox_p  = wilcox_pval,
      KS_p      = ks_pval,
      AD_p      = ad_pval,
      Loc       = L, 
      Scale     = S, 
      Shape     = G
    ))
    
    if(b == 28) {
      final_maxima <- block_maxima_temp
      final_fit <- gev_fit
      loc <- L; scale <- S; shape <- G
    }
  }, error = function(e) { message(paste("Error at block size", b, ":", e)) })
}

# Export Block Selection Table to LaTeX
results_latex <- results
colnames(results_latex) <- c("Block Size", "n Blks", "RMSE", "Runs (p)", "Wilcox (p)", "KS (p)", "AD (p)", "mu", "sigma", "xi")
print(xtable(results_latex, digits = 3), include.rownames = FALSE)

# ==========================================================
# 3. GEV BOOTSTRAP UNCERTAINTY (28-DAY BLOCKS)
# ==========================================================
boot_gev_lmom <- function(data, indices) {
  d <- data[indices]
  fit <- tryCatch(
    fevd(d, type = "GEV", method = "Lmoments"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(location = NA, scale = NA, shape = NA))
  p <- distill(fit)
  return(c(
    location = as.numeric(p["location"]),
    scale    = as.numeric(p["scale"]),
    shape    = as.numeric(p["shape"])
  ))
}

set.seed(123)
boot_gev_res <- boot(data = final_maxima, statistic = boot_gev_lmom, R = 2000)
gev_pars <- distill(final_fit)

gev_uncertainty_table <- data.frame(
  Parameter   = c("Location (μ)", "Scale (σ)", "Shape (ξ)"),
  Estimate    = round(as.numeric(gev_pars[c("location", "scale", "shape")]), 4),
  Std_Error   = round(apply(boot_gev_res$t, 2, sd, na.rm = TRUE), 4),
  Lower_95_CI = round(apply(boot_gev_res$t, 2, quantile, probs = 0.025, na.rm = TRUE), 4),
  Upper_95_CI = round(apply(boot_gev_res$t, 2, quantile, probs = 0.975, na.rm = TRUE), 4)
) %>%
  mutate(
    Formatted_Report = sprintf("%.4f (±%.4f)", Estimate, 1.96 * Std_Error),
    CI_95_Percent    = sprintf("[%.4f, %.4f]", Lower_95_CI, Upper_95_CI)
  )

cat("\n=== GEV L-MOMENTS PARAMETER UNCERTAINTY ===\n")
print(gev_uncertainty_table %>% select(Parameter, Estimate, Std_Error, CI_95_Percent), row.names = FALSE)

print(xtable(gev_uncertainty_table %>% select(Parameter, Estimate, Std_Error, CI_95_Percent),
             caption = "GEV parameter estimates ($\\\\mu, \\\\sigma, \\\\xi$) obtained via L-moments with 2,000 bootstrap standard errors and 95\\\\% confidence intervals.",
             label   = "tab:gev_lmom_uncertainty"),
      include.rownames = FALSE, booktabs = TRUE, comment = FALSE)

# ==========================================================
# 4. NON-STATIONARY GEV MODEL COMPARISON
# ==========================================================
n_blocks <- length(final_maxima)
time_cov <- 1:n_blocks
time_cov_sq <- time_cov^2

fit_stat <- fevd(final_maxima, type = "GEV")
fit_lin  <- fevd(final_maxima, type = "GEV", location.fun = ~time_cov, data = data.frame(time_cov))
fit_quad <- fevd(final_maxima, type = "GEV", location.fun = ~time_cov + time_cov_sq, data = data.frame(time_cov, time_cov_sq))

models <- list(fit_stat, fit_lin, fit_quad)
model_names <- c("Stationary", "Linear Trend", "Quadratic Trend")

comp_results <- data.frame(
  Model = model_names,
  Log_Lik = sapply(models, function(x) as.numeric(x$results$value) * -1),
  AIC = sapply(models, function(x) summary(x)$AIC),
  BIC = sapply(models, function(x) summary(x)$BIC),
  Par = c(3, 4, 5)
)

clean_comp <- as.data.frame(matrix(as.numeric(as.matrix(comp_results[,-1])), ncol = 4))
clean_comp <- cbind(Model = model_names, clean_comp)
colnames(clean_comp) <- c("Model", "Log-Lik.", "AIC", "BIC", "Par.")

print(xtable(clean_comp, 
             digits = c(0, 0, 2, 2, 2, 0), 
             caption = "Model Comparison for Non-Stationary GEV Analysis"), 
      include.rownames = FALSE, 
      comment = FALSE,
      booktabs = TRUE)

# ==========================================================
# 5. GEV PROBABILITY (P-P) PLOT & PRIMARY RETURN LEVELS
# ==========================================================
gev_par <- distill(final_fit)
mu_hat    <- as.numeric(gev_par["location"])
sigma_hat <- as.numeric(gev_par["scale"])
xi_hat    <- as.numeric(gev_par["shape"])

obs <- sort(final_maxima)
empirical <- ppoints(length(obs))
theoretical <- pgev(obs, loc = mu_hat, scale = sigma_hat, shape = xi_hat)

pp_data <- data.frame(Theoretical = theoretical, Empirical = empirical)

gev_pp_plot <- ggplot(pp_data, aes(x = Theoretical, y = Empirical)) +
  geom_point(colour = "steelblue", size = 3) +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed", linewidth = 0.8) +
  labs(title = "Probability Plot for the Fitted GEV Model", x = "Theoretical Probability", y = "Empirical Probability") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))

print(gev_pp_plot)

# Return Level Grid & Figures (2 to 20 Years)
blocks_per_year <- 365.25 / 28
n_obs <- length(obs)
empirical_ranks <- 1:n_obs
empirical_periods_blocks <- (n_obs + 0.12) / (n_obs + 0.38 - empirical_ranks)
empirical_periods_years  <- empirical_periods_blocks / blocks_per_year

empirical_df <- data.frame(Years = empirical_periods_years, Observed = obs) %>% filter(Years <= 20)
min_year <- min(empirical_df$Years, na.rm = TRUE)

primary_years <- c(2, 5, 10, 20)
grid_periods_years <- exp(seq(log(min_year), log(20), length.out = 200))

boot_gev_combined <- function(data, indices, grid_years, discrete_years, B) {
  d <- data[indices]
  fit <- tryCatch(fevd(d, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit)) return(c(rep(NA, length(grid_years)), rep(NA, length(discrete_years))))
  rl_grid <- as.numeric(return.level(fit, return.period = grid_years * B))
  rl_disc <- as.numeric(return.level(fit, return.period = discrete_years * B))
  return(c(rl_grid, rl_disc))
}

set.seed(123)
boot_combined_res <- boot(
  data = final_maxima, 
  statistic = boot_gev_combined, 
  R = 1000, 
  grid_years = grid_periods_years,
  discrete_years = primary_years,
  B = blocks_per_year
)

n_grid <- length(grid_periods_years)
smooth_plot_df <- data.frame(
  Years        = grid_periods_years,
  Estimated_RL = boot_combined_res$t0[1:n_grid],
  Lower_CI     = pmax(0, apply(boot_combined_res$t[, 1:n_grid], 2, quantile, probs = 0.025, na.rm = TRUE)),
  Upper_CI     = apply(boot_combined_res$t[, 1:n_grid], 2, quantile, probs = 0.975, na.rm = TRUE)
)

n_disc <- length(primary_years)
disc_indices <- (n_grid + 1):(n_grid + n_disc)

primary_summary_plot <- data.frame(
  Years        = primary_years,
  Estimated_RL = boot_combined_res$t0[disc_indices],
  Lower_CI     = pmax(0, apply(boot_combined_res$t[, disc_indices], 2, quantile, probs = 0.025, na.rm = TRUE)),
  Upper_CI     = apply(boot_combined_res$t[, disc_indices], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# Construction of `gev_return_table` for Section 9 LaTeX Output
gev_return_table <- primary_summary_plot %>%
  transmute(
    `Return Period` = paste0(Years, "-year"),
    Estimate        = round(Estimated_RL, 2),
    `Lower 95% CI`  = round(Lower_CI, 2),
    `Upper 95% CI`  = round(Upper_CI, 2)
  )

main_gev_rl_plot <- ggplot() +
  geom_ribbon(data = smooth_plot_df, aes(x = Years, ymin = Lower_CI, ymax = Upper_CI), fill = "steelblue", alpha = 0.2) +
  geom_line(data = smooth_plot_df, aes(x = Years, y = Lower_CI), color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = smooth_plot_df, aes(x = Years, y = Upper_CI), color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = smooth_plot_df, aes(x = Years, y = Estimated_RL), color = "darkred", linewidth = 0.9) +
  geom_point(data = primary_summary_plot, aes(x = Years, y = Estimated_RL), color = "darkred", size = 2.5) +
  geom_point(data = empirical_df, aes(x = Years, y = Observed), shape = 21, fill = "white", color = "black", size = 2, stroke = 0.8) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20), labels = c("0.5", "1", "2", "5", "10", "20")) +
  coord_cartesian(xlim = c(min_year * 0.9, 22), ylim = c(0, NA)) +
  labs(x = "Return Period (Years)", y = "Return Level (Cases per Million)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = "plain"))

print(main_gev_rl_plot)

# ==========================================================
# 6. POT DIAGNOSTICS & THRESHOLD SELECTION
# ==========================================================
active_cases <- zambia_data %>% filter(!is.na(cases), cases > 0) %>% pull(cases)
n_positive_obs <- length(active_cases)
study_days <- as.numeric(max(zambia_data$Day, na.rm = TRUE) - min(zambia_data$Day, na.rm = TRUE)) + 1
study_years <- study_days / 365.25

percentiles <- seq(0.85, 0.95, by = 0.01)
u_candidates <- quantile(active_cases, probs = percentiles, na.rm = TRUE)

evaluate_threshold <- function(u, percentile_label, data, n_total, years, run_length = 7) {
  exceedances <- data[data > u]
  raw_N       <- length(exceedances)
  if (raw_N == 0) return(NULL)
  
  declustered_series <- decluster(x = data, threshold = u, method = "runs", r = run_length)
  cluster_N <- sum(declustered_series > u, na.rm = TRUE)
  theta <- cluster_N / raw_N
  raw_rate_daily  <- raw_N / n_total
  lambda_c_annual <- cluster_N / years
  
  mean_excess <- mean(exceedances - u)
  se_excess   <- sd(exceedances - u) / sqrt(raw_N)
  
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
    Threshold_u          = round(u, 2),
    Raw_N                = raw_N,
    Clusters_Nc          = cluster_N,
    Theta                = round(theta, 3),
    Raw_Rate_Daily       = round(raw_rate_daily, 4),
    Cluster_Rate_Annual  = round(lambda_c_annual, 4),
    Mean_Excess          = mean_excess,
    Lower_CI             = pmax(0, mean_excess - 1.96 * se_excess),
    Upper_CI             = mean_excess + 1.96 * se_excess,
    Scale_sigma          = round(sigma_hat, 4),
    Shape_xi             = round(xi_hat, 4)
  )
}

threshold_table_full <- map2_df(u_candidates, paste0(round(percentiles * 100), "%"), 
                                ~ evaluate_threshold(.x, .y, active_cases, n_positive_obs, study_years, 7))

# Diagnostic Plots
mrl_gg <- ggplot(threshold_table_full, aes(x = Threshold_u, y = Mean_Excess)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(aes(size = Clusters_Nc), color = "darkred") +
  geom_vline(xintercept = 61.26, linetype = "dashed") +
  labs(x = "Threshold (u)", y = "Mean Excess") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

shape_stab_gg <- ggplot(threshold_table_full, aes(x = Threshold_u, y = Shape_xi)) +
  geom_ribbon(aes(ymin = Shape_xi - 0.1, ymax = Shape_xi + 0.1), fill = "darkgreen", alpha = 0.1) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2) +
  geom_vline(xintercept = 61.26, linetype = "dashed", color = "red") +
  labs(x = "Threshold (u)", y = "Shape Estimate (ξ)") +
  theme_bw() + theme(panel.grid = element_blank())

scale_stab_gg <- ggplot(threshold_table_full, aes(x = Threshold_u, y = Scale_sigma)) +
  geom_ribbon(aes(ymin = Scale_sigma * 0.9, ymax = Scale_sigma * 1.1), fill = "darkblue", alpha = 0.1) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 2) +
  geom_vline(xintercept = 61.26, linetype = "dashed", color = "red") +
  labs(x = "Threshold (u)", y = "Scale Estimate (σ)") +
  theme_bw() + theme(panel.grid = element_blank())

print(mrl_gg)
print(shape_stab_gg)
print(scale_stab_gg)

# ==========================================================
# 7. SELECTED POT CLUSTER & SENSITIVITY ANALYSIS
# ==========================================================
u_selected <- 61.26
run_length <- 7

exceedance_data <- zambia_data %>%
  filter(cases > u_selected) %>%
  arrange(Day) %>%
  mutate(
    gap_days    = as.numeric(Day - lag(Day)),
    new_cluster = ifelse(is.na(gap_days) | gap_days > run_length, 1, 0),
    cluster_id  = cumsum(new_cluster)
  )

cluster_data <- exceedance_data %>%
  group_by(cluster_id) %>%
  summarise(
    Cluster_Start   = min(Day),
    Cluster_End     = max(Day),
    Duration_Days   = as.numeric(Cluster_End - Cluster_Start) + 1,
    Cluster_Size    = n(),
    Cluster_Maximum = max(cases),
    .groups         = "drop"
  )

raw_N_selected        <- nrow(exceedance_data)
cluster_N_selected    <- nrow(cluster_data)
theta_selected        <- cluster_N_selected / raw_N_selected
mean_cluster_duration <- mean(cluster_data$Duration_Days)
mean_cluster_size     <- mean(cluster_data$Cluster_Size)
lambda_cluster_annual <- cluster_N_selected / study_years

cluster_peaks <- cluster_data$Cluster_Maximum
final_gpd_fit <- fevd(x = cluster_peaks, threshold = u_selected, type = "GP", method = "Lmoments")

gpd_parameters <- distill(final_gpd_fit)
sigma_hat      <- as.numeric(gpd_parameters["scale"])
xi_hat          <- as.numeric(gpd_parameters["shape"])

run_lengths <- c(1, 2, 3, 4, 5, 7, 10, 14, 21)
run_length_sensitivity <- map_df(run_lengths, function(r) {
  dec <- decluster(x = active_cases, threshold = u_selected, method = "runs", r = r)
  cluster_N <- sum(dec > u_selected, na.rm = TRUE)
  theta     <- cluster_N / raw_N_selected
  lambda_c  <- cluster_N / study_years
  data.frame(
    Run_Length_Days      = r,
    Raw_Exceedances      = raw_N_selected,
    Independent_Clusters = cluster_N,
    Theta                = round(theta, 4),
    Annual_Cluster_Rate  = round(lambda_c, 4)
  )
})

# POT Bootstrap Uncertainty
set.seed(123)
cluster_sizes <- cluster_data$Cluster_Size
theta_boot <- replicate(2000, {
  sampled_sizes <- sample(cluster_sizes, size = length(cluster_sizes), replace = TRUE)
  length(sampled_sizes) / sum(sampled_sizes)
})
theta_CI <- quantile(theta_boot, probs = c(0.025, 0.975), na.rm = TRUE)

gpd_boot <- map_dfr(1:2000, function(i) {
  sampled_peaks <- sample(cluster_peaks, size = length(cluster_peaks), replace = TRUE)
  fit_boot <- tryCatch(fevd(x = sampled_peaks, threshold = u_selected, type = "GP", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit_boot)) return(data.frame(Sigma = NA_real_, Xi = NA_real_))
  pars_boot <- distill(fit_boot)
  data.frame(Sigma = as.numeric(pars_boot["scale"]), Xi = as.numeric(pars_boot["shape"]))
})

sigma_CI <- quantile(gpd_boot$Sigma[is.finite(gpd_boot$Sigma)], probs = c(0.025, 0.975), na.rm = TRUE)
xi_CI    <- quantile(gpd_boot$Xi[is.finite(gpd_boot$Xi)], probs = c(0.025, 0.975), na.rm = TRUE)

parameter_uncertainty <- data.frame(
  Parameter   = c("Theta (θ)", "Scale (σ)", "Shape (ξ)"),
  Estimate    = round(c(theta_selected, sigma_hat, xi_hat), 4),
  Lower_95_CI = round(c(theta_CI[1], sigma_CI[1], xi_CI[1]), 4),
  Upper_95_CI = round(c(theta_CI[2], sigma_CI[2], xi_CI[2]), 4)
)

# Time Series Exceedances Plot
plot_ts_df <- zambia_data %>% mutate(Status = ifelse(cases > u_selected, "Exceedance", "Baseline"))
ts_gg <- ggplot(plot_ts_df, aes(x = Day, y = cases)) +
  geom_line(color = "gray80", linewidth = 0.5) +
  geom_point(data = filter(plot_ts_df, Status == "Exceedance"), aes(color = Status), size = 1.5) +
  geom_hline(yintercept = u_selected, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("Exceedance" = "darkred")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", limits = c(as.Date("2020-01-01"), as.Date("2023-12-31"))) +
  labs(title = "COVID-19 Exceedances in Zambia", x = "Timeline (2020 - 2023)", y = "Active Cases (per Million)") +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank())

print(ts_gg)

# GPD P-P Plot
exceedances_gpd <- sort(cluster_peaks - u_selected)
n_exc <- length(exceedances_gpd)
empirical_gpd <- ppoints(n_exc)
theoretical_gpd <- pgpd(exceedances_gpd, loc = 0, scale = sigma_hat, shape = xi_hat)

pp_plot <- ggplot(data.frame(Empirical = empirical_gpd, Theoretical = theoretical_gpd), aes(x = Theoretical, y = Empirical)) +
  geom_point(size = 3, colour = "steelblue") +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed", linewidth = 0.8) +
  labs(title = "Probability Plot for the Fitted GPD Model", x = "Theoretical Probability", y = "Empirical Probability") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))

print(pp_plot)

# ==========================================================
# 8. GPD RETURN LEVEL BOOTSTRAP (95% CI)
# ==========================================================
return_periods <- c(2, 5, 10, 20)
lambda_c <- lambda_cluster_annual

calc_gpd_rl <- function(T_val, u, sigma, xi, rate) {
  if (abs(xi) < 1e-8) {
    return(u + sigma * log(rate * T_val))
  } else {
    return(u + (sigma / xi) * ((rate * T_val)^xi - 1))
  }
}

estimated_RL <- sapply(return_periods, calc_gpd_rl, u = u_selected, sigma = sigma_hat, xi = xi_hat, rate = lambda_c)

gpd_boot_stat <- function(data, indices) {
  sample_peaks <- data[indices]
  if (length(unique(sample_peaks)) < 3) return(rep(NA_real_, length(return_periods)))
  fit_b <- tryCatch(fevd(x = sample_peaks, threshold = u_selected, type = "GP", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit_b)) return(rep(NA_real_, length(return_periods)))
  pars_b <- tryCatch(distill(fit_b), error = function(e) NULL)
  if (is.null(pars_b)) return(rep(NA_real_, length(return_periods)))
  sigma_b <- as.numeric(pars_b["scale"]); xi_b <- as.numeric(pars_b["shape"])
  if (is.na(sigma_b) || is.na(xi_b) || sigma_b <= 0) return(rep(NA_real_, length(return_periods)))
  sapply(return_periods, calc_gpd_rl, u = u_selected, sigma = sigma_b, xi = xi_b, rate = lambda_c)
}

set.seed(123)
boot_results <- boot(data = cluster_peaks, statistic = gpd_boot_stat, R = 5000)

valid_rows <- complete.cases(boot_results$t)
valid_boot <- boot_results$t[valid_rows, , drop = FALSE]

cat(
  "Successful Bootstrap Replicates:",
  sum(valid_rows),
  "out of",
  nrow(boot_results$t),
  "\n\n"
)

lower_CI <- apply(valid_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
upper_CI <- apply(valid_boot, 2, quantile, probs = 0.975, na.rm = TRUE)

return_level_table <- data.frame(
  Return_Period = paste0(return_periods, "-year"),
  Lower_95_CI   = round(lower_CI, 2),
  Estimate      = round(estimated_RL, 2),
  Upper_95_CI   = round(upper_CI, 2)
)

kable(
  return_level_table,
  caption = "POT-GPD Annual Return Level Estimates with 95% Bootstrap CIs",
  col.names = c("Return Period", "Lower 95% CI", "Estimate", "Upper 95% CI"),
  align = c("l", "r", "r", "r")
)

latex_rl_table <- xtable(
  return_level_table,
  caption = "POT-GPD Return level estimates (cases per million) with 95\\% confidence intervals.",
  label   = "tab:gpd_return_levels",
  digits  = c(0, 0, 2, 2, 2)
)

print(
  latex_rl_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  comment = FALSE
)

# Construction of `long_return_table` for Sensitivity Analysis (50 & 100 Years)
long_years <- c(50, 100)
boot_gev_long <- function(data, indices, periods, B) {
  d <- data[indices]
  fit <- tryCatch(fevd(d, type = "GEV", method = "Lmoments"), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, length(periods)))
  return(as.numeric(return.level(fit, return.period = periods * B)))
}

set.seed(123)
boot_long_res <- boot(data = final_maxima, statistic = boot_gev_long, R = 1000, periods = long_years, B = blocks_per_year)

long_return_table <- data.frame(
  `Return Period` = paste0(long_years, "-year"),
  Estimate        = round(boot_long_res$t0, 2),
  `Lower 95% CI`  = round(pmax(0, apply(boot_long_res$t, 2, quantile, probs = 0.025, na.rm = TRUE)), 2),
  `Upper 95% CI`  = round(apply(boot_long_res$t, 2, quantile, probs = 0.975, na.rm = TRUE), 2)
)

# ==========================================================
# 9. PUBLICATION-READY LATEX TABLES
# ==========================================================
cat("\n==========================================================\n")
cat("LATEX: CLUSTER TABLE\n")
cat("==========================================================\n")
print(xtable(cluster_data, caption = "Identified exceedance clusters at the selected threshold using a 7-day runs parameter.", label = "tab:cluster_details"), include.rownames = FALSE, booktabs = TRUE, comment = FALSE, digits = 2)

cat("\n==========================================================\n")
cat("LATEX: RUN-LENGTH SENSITIVITY TABLE\n")
cat("==========================================================\n")
print(xtable(run_length_sensitivity, caption = "Sensitivity of runs declustering to alternative run lengths.", label = "tab:run_length_sensitivity"), include.rownames = FALSE, booktabs = TRUE, comment = FALSE, digits = 4)

cat("\n==========================================================\n")
cat("LATEX: PARAMETER UNCERTAINTY TABLE\n")
cat("==========================================================\n")
print(xtable(parameter_uncertainty, caption = "Bootstrap uncertainty estimates for the extremal index and GPD parameters at the selected threshold.", label = "tab:pot_uncertainty"), include.rownames = FALSE, booktabs = TRUE, comment = FALSE, digits = 4)

cat("\n==========================================================\n")
cat("LATEX: PRIMARY GEV RETURN LEVEL TABLE\n")
cat("==========================================================\n")
print(xtable(gev_return_table, caption = "Primary stationary GEV return levels for 2-, 5-, 10- and 20-year return periods.", label = "tab:gev_return_levels"), include.rownames = FALSE, booktabs = TRUE, comment = FALSE, digits = 2)

cat("\n==========================================================\n")
cat("LATEX: LONG-HORIZON SENSITIVITY TABLE\n")
cat("==========================================================\n")
print(xtable(long_return_table, caption = "Long-horizon stationary GEV return levels presented as sensitivity analysis.", label = "tab:long_horizon_return_levels"), include.rownames = FALSE, booktabs = TRUE, comment = FALSE, digits = 2)

cat("\n==========================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("==========================================================\n")
