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








#====================================================================
# 1. Keep only positive cases
#--------------------------------------------------------------------
active_cases <- cases_vector[!is.na(cases_vector) & cases_vector > 0]

#-------------------------------------------------
# 2. Candidate thresholds
#-------------------------------------------------
percentiles <- seq(0.85, 0.95, by = 0.01)
u_values <- quantile(active_cases, percentiles)

#-------------------------------------------------
# 3. Total number of observations
#-------------------------------------------------
n_total <- length(active_cases)

#-------------------------------------------------
# 4. Create threshold comparison table
#-------------------------------------------------

threshold_table <- map_df(seq_along(u_values), function(i){
  
  u <- as.numeric(u_values[i])
  
  ## Raw exceedances
  exceed <- active_cases[active_cases > u]
  raw_N <- length(exceed)
  
  ## Declustering (7-day run)
  dec <- decluster(active_cases,
                   threshold = u,
                   method = "runs",
                   r = 7)
  
  ## Independent clusters
  cluster_N <- sum(dec > u, na.rm = TRUE)
  
  ## Extremal Index
  theta <- ifelse(raw_N > 0,
                  cluster_N/raw_N,
                  NA)
  
  ## Rate
  lambda <- cluster_N/n_total
  
  ## Fit GPD
  fit <- try(
    fevd(dec,
         threshold = u,
         type = "GP",
         method = "Lmoments"),
    silent = TRUE
  )
  
  if(inherits(fit,"try-error")){
    
    sigma <- NA
    xi <- NA
    
  } else{
    
    pars <- distill(fit)
    
    sigma <- as.numeric(pars["scale"])
    xi <- as.numeric(pars["shape"])
    
  }
  
  data.frame(
    
    Percentile = paste0(round(percentiles[i]*100),"%"),
    
    Threshold_u = round(u,2),
    
    Raw_N = raw_N,
    
    Clusters_Nc = cluster_N,
    
    Theta = round(theta,3),
    
    Rate_lambda = round(lambda,4),
    
    Scale_sigma = round(sigma,4),
    
    Shape_xi = round(xi,4)
    
  )
  
})

#-------------------------------------------------
# Display
#-------------------------------------------------

print(threshold_table,row.names=FALSE)

#-----------------------------------------------
# Return Level Table
#-----------------------------------------------

# Return periods (years)
return_periods <- c(2, 5, 10, 20, 50, 100)

# Compute return levels with confidence intervals
rl <- ci(
  fit_gpd,
  type = "return.level",
  return.period = return_periods,
  method = "normal"
)

# Create publication-quality table
return_level_table <- data.frame(
  `Return Period` = paste0(return_periods, "-year"),
  `Lower 95% CI`  = round(rl[,1], 2),
  Estimate        = round(rl[,2], 2),
  `Upper 95% CI`  = round(rl[,3], 2)
)

# Display table
print(return_level_table, row.names = FALSE)








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













