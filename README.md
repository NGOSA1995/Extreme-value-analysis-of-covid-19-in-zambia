
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
    gev_fit <- fevd(block_maxima_temp, type = "GEV", method = "Lmoments")
    params <- distill(gev_fit)
    
    # Validation Stats
    p_seq <- ppoints(length(block_maxima_temp))
    q_theory <- qgev(p_seq, loc = params["location"], scale = params["scale"], shape = params["shape"])
    rmse_val <- sqrt(mean((sort(block_maxima_temp) - q_theory)^2))
    ks_pval <- ks.test(block_maxima_temp, "pgev", loc = params["location"], 
                       scale = params["scale"], shape = params["shape"])$p.value
    
    results <- rbind(results, data.frame(
      BlockSize = b, nBlocks = length(block_maxima_temp), RMSE = rmse_val, KS_pval = ks_pval,
      Location = params["location"], Scale = params["scale"], Shape = params["shape"]
    ))
    
    # Store 28-day variables for final plotting
    if(b == 28) {
      final_maxima <- block_maxima_temp
      final_fit <- gev_fit
      loc <- params["location"]; scale <- params["scale"]; shape <- params["shape"]
    }
  }, error = function(e) {})
}

print(xtable(results, digits=3), include.rownames = FALSE)

# ==========================================================
#  FIGURE: Q-Q PLOT 
# ==========================================================
obs_sorted <- sort(final_maxima)
theory_quantiles <- qgev(ppoints(length(final_maxima)), loc = loc, scale = scale, shape = shape)
qq_data <- data.frame(Observed = obs_sorted, Theoretical = theory_quantiles)

ggplot(qq_data, aes(x = Theoretical, y = Observed)) +
  geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 0.7) + 
  geom_point(shape = 21, fill = "steelblue", color = "white", size = 3, alpha = 0.8) +
  labs(x = "Theoretical Quantiles", y = "Observed Maxima (Cases/Million)",
       subtitle = NULL) +
  theme_bw(base_size = 12) + 
  theme(
    panel.grid.major = element_blank(), # Grid removed
    panel.grid.minor = element_blank(), # Grid removed
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "italic")
  )

# ==========================================================
# BOOTSTRAP RETURN LEVELS
# ==========================================================
get_return_levels <- function(data, indices, periods) {
  d <- data[indices]
  fit <- tryCatch({ fevd(d, type = "GEV", method = "Lmoments") }, error = function(e) NULL)
  if(is.null(fit)) return(rep(NA, length(periods)))
  return(as.numeric(return.level(fit, return.period = periods)))
}

return_periods <- c(2, 5, 10, 20, 50, 100)
set.seed(123)
boot_results <- boot(data = final_maxima, statistic = get_return_levels, R = 1000, periods = return_periods)

boot_summary <- data.frame(
  Return_Period = return_periods,
  Estimated_RL  = boot_results$t0,
  Lower_CI      = apply(boot_results$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  Upper_CI      = apply(boot_results$t, 2, quantile, probs = 0.975, na.rm = TRUE)
)

# ==========================================================
#  FIGURE: RETURN LEVEL PLOT 
# ==========================================================
ggplot(boot_summary, aes(x = Return_Period, y = Estimated_RL)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "gray80", alpha = 0.5) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(color = "black", size = 2, shape = 16) +
  scale_x_log10(breaks = return_periods) +
  labs(x = "Return Period (Years)", y = "Return Level (Cases per Million)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(), # Grid removed
    panel.grid.minor = element_blank(), # Grid removed
    axis.title = element_text(face = "plain")
  )


# ==========================================================
# COMPLETE EXTREME VALUE ANALYSIS: ZAMBIA COVID-19
# ==========================================================
library(extRemes)
library(tidyverse)

# ----------------------------------------------------------
#  DATA ALIGNMENT
# ----------------------------------------------------------
if (!exists("dates")) {
  dates <- seq(as.Date("2020-03-18"), by = "day", length.out = length(active_cases))
}

u_selected <- 61.26

# ----------------------------------------------------------
#  THRESHOLD SENSITIVITY DATA GENERATION
# ----------------------------------------------------------
percentiles <- seq(0.85, 0.95, by = 0.01)
u_values <- quantile(active_cases, percentiles)
n_total <- length(active_cases)

threshold_table_full <- map_df(seq_along(u_values), function(i) {
  u <- as.numeric(u_values[i])
  exc <- active_cases[active_cases > u]
  n_exc <- length(exc)
  
  m_excess <- mean(exc - u, na.rm = TRUE)
  s_error  <- sd(exc - u, na.rm = TRUE) / sqrt(n_exc)
  
  data_dec_tmp <- decluster(active_cases, threshold = u, method = "runs", r = 7)
  n_clusters <- sum(data_dec_tmp > u, na.rm = TRUE)
  
  fit_tmp <- try(fevd(data_dec_tmp, threshold = u, type = "GP", method = "Lmoments"), silent = TRUE)
  
  if(!inherits(fit_tmp, "try-error")) {
    params <- distill(fit_tmp)
    data.frame(
      Percentile = names(u_values)[i],
      Threshold_u = u,
      Mean_Excess = m_excess,
      Lower_CI = m_excess - (1.96 * s_error),
      Upper_CI = m_excess + (1.96 * s_error),
      Clusters_Nc = n_clusters,
      Scale_sigma = as.numeric(params["scale"]),
      Shape_xi = as.numeric(params["shape"])
    )
  }
})

# ----------------------------------------------------------
#  FINAL MODEL FIT (L-MOMENTS)
# ----------------------------------------------------------
declustered_data <- decluster(active_cases, threshold = u_selected, method = "runs", r = 7)

final_gpd_fit <- fevd(declustered_data, 
                      threshold = u_selected, 
                      type = "GP", 
                      method = "Lmoments")

# ----------------------------------------------------------
#  ANNUAL RETURN LEVELS (TABLE)
# ----------------------------------------------------------
years <- c(1, 2, 5, 10, 25, 50, 100)
long_term_rl <- return.level(final_gpd_fit, return.period = years * 365.25, do.ci = TRUE)

annual_table <- data.frame(
  Return_Period_Years = years,
  Lower_95_CI = round(as.numeric(long_term_rl[,1]), 2),
  Estimate = round(as.numeric(long_term_rl[,2]), 2),
  Upper_95_CI = round(as.numeric(long_term_rl[,3]), 2)
)

print("--- ANNUAL RETURN LEVELS (1 TO 100 YEARS) ---")
print(annual_table)

# ----------------------------------------------------------
#  CONSOLIDATED FIGURES (NO GRIDS)
# ----------------------------------------------------------

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
# NON-STATIONARY GEV MODEL COMPARISON
# ==========================================================

#  Prepare data and covariate (using 28-day maxima as example)
# 'final_maxima' should be your vector of block maxima
n_blocks <- length(final_maxima)
time_cov <- 1:n_blocks
time_cov_sq <- time_cov^2

#  Fit the three candidate models
# Stationary Model
fit_stat <- fevd(final_maxima, type = "GEV")

# Linear Trend Model (Location parameter depends on Time)
fit_lin  <- fevd(final_maxima, type = "GEV", 
                 location.fun = ~time_cov, data = data.frame(time_cov))

# Quadratic Trend Model (Location parameter depends on Time + Time^2)
fit_quad <- fevd(final_maxima, type = "GEV", 
                 location.fun = ~time_cov + time_cov_sq, 
                 data = data.frame(time_cov, time_cov_sq))

#  Extract Metrics and Construct Table
models <- list(fit_stat, fit_lin, fit_quad)
model_names <- c("Stationary", "Linear Trend", "Quadratic Trend")

comp_results <- data.frame(
  Model = model_names,
  Log_Lik = sapply(models, function(x) as.numeric(x$results$value) * -1), # Convert to positive log-lik
  AIC = sapply(models, function(x) summary(x)$AIC),
  BIC = sapply(models, function(x) summary(x)$BIC),
  Par = c(3, 4, 5) # Number of parameters for each model
)

#  Export to LaTeX using xtable
library(xtable)

# Ensure numeric cleanliness
clean_comp <- as.data.frame(matrix(as.numeric(as.matrix(comp_results[,-1])), ncol = 4))
clean_comp <- cbind(Model = model_names, clean_comp)
colnames(clean_comp) <- c("Model", "Log-Lik.", "AIC", "BIC", "Par.")

print(xtable(clean_comp, 
             digits = c(0, 0, 2, 2, 2, 0), 
             caption = "Model Comparison for Non-Stationary GEV Analysis"), 
      include.rownames = FALSE, 
      comment = FALSE,
      booktabs = TRUE)
