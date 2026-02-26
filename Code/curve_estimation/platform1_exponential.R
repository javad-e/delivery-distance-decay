library(tidyverse)
library(lme4)
library(splines)
library(scales)

CONTROL_VARIABLES <- c(
  "num_vendors_dest",
  "income_cat_origin",
  "income_cat_dest",
  "vendor_entropy_origin",
  "vendor_entropy_destination",
  "pct_incentivized",
  "population_dest",
  "avg_delivery_fee"
)


run_spline_model_return_curve <- function(df, 
                                          control_vars = CONTROL_VARIABLES,
                                          dist_min = 0, 
                                          dist_max = 18, 
                                          spline_step = 2, 
                                          pred_points = 200) {
  # Filter and prepare data
  df <- df %>%
    filter(euclidean_km >= dist_min, euclidean_km <= dist_max) %>%
    mutate(
      origin = factor(origin),
      destination = factor(destination),
      log_flow = log(flow + 1)
    )
  
  # Check which control variables are actually present in the data
  available_controls <- intersect(control_vars, names(df))
  missing_controls <- setdiff(control_vars, names(df))
  
  if (length(missing_controls) > 0) {
    warning("The following control variables are not in the dataset and will be skipped: ",
            paste(missing_controls, collapse = ", "))
  }
  
  if (length(available_controls) > 0) {
    cat("Using control variables:", paste(available_controls, collapse = ", "), "\n")
  } else {
    cat("No control variables will be used\n")
  }
  
  # Standardize control variables (for better convergence)
  for (var in available_controls) {
    if (is.numeric(df[[var]])) {
      df[[paste0(var, "_std")]] <- scale(df[[var]])[, 1]
    } else {
      warning(paste("Control variable", var, "is not numeric and will be skipped"))
      available_controls <- setdiff(available_controls, var)
    }
  }
  
  # Create spline basis
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- ns(
    df$euclidean_km,
    knots = knots,
    intercept = FALSE,
    Boundary.knots = c(dist_min, dist_max)
  )
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  
  df <- bind_cols(df, as.data.frame(spline_basis))
  
  # Build formula with spline terms and control variables
  spline_terms <- paste(paste0("spline_", seq_len(ncol(spline_basis)) - 1), 
                        collapse = " + ")
  
  control_terms <- ""
  if (length(available_controls) > 0) {
    control_terms_std <- paste0(available_controls, "_std")
    control_terms <- paste0(" + ", paste(control_terms_std, collapse = " + "))
  }
  
  formula_text <- paste0("log_flow ~ 1 + ", spline_terms, control_terms,
                         " + (1 | origin) + (1 | destination)")
  
  cat("Model formula:", formula_text, "\n\n")
  
  # Fit mixed-effects model
  model <- lmer(as.formula(formula_text), data = df)
  
  # Generate predictions on grid (partial effect: fixed effects only)
  distance_grid <- tibble(distance = seq(dist_min, dist_max, length.out = pred_points))
  spline_basis_grid <- ns(
    distance_grid$distance,
    knots = knots,
    intercept = FALSE,
    Boundary.knots = c(dist_min, dist_max)
  )
  colnames(spline_basis_grid) <- paste0("spline_", seq_len(ncol(spline_basis_grid)) - 1)
  
  distance_grid <- bind_cols(distance_grid, as.data.frame(spline_basis_grid))
  distance_grid$origin <- factor(levels(df$origin)[1], levels = levels(df$origin))
  distance_grid$destination <- factor(levels(df$destination)[1], levels = levels(df$destination))
  
  for (var in available_controls) {
    distance_grid[[paste0(var, "_std")]] <- 0
  }
  
  X_grid <- model.matrix(
    reformulate(c(paste0("spline_", seq_len(ncol(spline_basis)) - 1), 
                  paste0(available_controls, "_std"))),
    data = distance_grid
  )
  
  # Get fixed effects coefficients and variance-covariance matrix
  beta <- fixef(model)
  vcov_beta <- as.matrix(vcov(model))
  
  # Calculate predictions
  distance_grid$pred_log_flow <- as.vector(X_grid %*% beta)
  
  # Calculate standard errors from parameter uncertainty only
  distance_grid$se_log_flow <- sqrt(diag(X_grid %*% vcov_beta %*% t(X_grid)))
  
  # Calculate confidence intervals on log scale, then transform
  distance_grid <- distance_grid %>%
    mutate(
      log_ci_lower = pred_log_flow - 1.96 * se_log_flow,
      log_ci_upper = pred_log_flow + 1.96 * se_log_flow,
      pred_flow = pmax(exp(pred_log_flow) - 1, 0),
      ci_lower = pmax(exp(log_ci_lower) - 1, 0),
      ci_upper = pmax(exp(log_ci_upper) - 1, 0)
    )
  
  # Calculate R² (conditional)
  df$pred_log_flow_full <- predict(model, re.form = NULL)
  ss_total <- sum((df$log_flow - mean(df$log_flow))^2)
  ss_res <- sum((df$log_flow - df$pred_log_flow_full)^2)
  r2_conditional <- 1 - ss_res / ss_total
  
  # Return partial effect curve and raw data with parameter uncertainty CIs
  list(
    partial_effect = distance_grid %>%
      dplyr::select(distance, pred_flow, ci_lower, ci_upper) %>%
      mutate(r2 = r2_conditional),
    raw_data = df %>%
      dplyr::select(distance = euclidean_km, flow),
    model = model,
    control_vars_used = available_controls
  )
}


# Fit exponential decay model to partial effect predictions
fit_exponential_model <- function(partial_effect_data, dist_min = 0, dist_max = 18, pred_points = 200) {
  # Fit exponential decay to the spline partial effect predictions
  # flow = a * exp(-b * distance)
  # Using nonlinear least squares
  tryCatch({
    model_exp <- nls(
      pred_flow ~ a * exp(-b * distance),
      data = partial_effect_data,
      start = list(a = max(partial_effect_data$pred_flow), b = 0.1),
      control = nls.control(maxiter = 1000)
    )
    
    # Generate predictions
    distance_grid <- tibble(distance = seq(dist_min, dist_max, length.out = pred_points))
    distance_grid$pred_flow_exp <- predict(model_exp, newdata = distance_grid)
    
    # Calculate R² (fit to partial effect curve)
    ss_total <- sum((partial_effect_data$pred_flow - mean(partial_effect_data$pred_flow))^2)
    ss_res <- sum((partial_effect_data$pred_flow - predict(model_exp, newdata = partial_effect_data))^2)
    r2_exp <- 1 - ss_res / ss_total
    
    cat("\n=== EXPONENTIAL MODEL ===\n")
    cat("Formula: flow = a * exp(-b * distance)\n")
    cat("Parameters:\n")
    print(coef(model_exp))
    cat("R²:", r2_exp, "\n\n")
    
    list(
      predictions = distance_grid %>%
        mutate(r2 = r2_exp),
      model = model_exp,
      r2 = r2_exp
    )
  }, error = function(e) {
    cat("Error fitting exponential model:", e$message, "\n")
    return(NULL)
  })
}


# --- Process Dubai Data ---
cat("\n=== Processing Dubai Dataset ===\n")
df <- read_csv("~/imperial/Reach of Last Mile/platform1_dubai_curve_estimation.csv", show_col_types = FALSE)
cat("Dataset dimensions:", nrow(df), "rows,", ncol(df), "columns\n")

results <- run_spline_model_return_curve(df, 
                                         control_vars = CONTROL_VARIABLES,
                                         dist_max = 18, 
                                         spline_step = 2)

# Extract results
express_partial <- results$partial_effect
express_raw <- results$raw_data

# Fit exponential model to the partial effect curve
exp_results <- fit_exponential_model(express_partial, dist_max = 18)

# --- Display Model Summary ---
cat("\n\n=== MODEL SUMMARY ===\n")
cat("Control variables used:", paste(results$control_vars_used, collapse = ", "), "\n")
cat("Spline Model R² (Conditional):", unique(express_partial$r2), "\n")
if (!is.null(exp_results)) {
  cat("Exponential Model R²:", exp_results$r2, "\n")
}

# --- Plot 1: Partial Effect Curve with 95% CI (Original Scale) with Exponential ---
if (!is.null(exp_results)) {
  p1 <- ggplot() +
    # Spline confidence interval
    geom_ribbon(
      data = express_partial,
      aes(x = distance, ymin = ci_lower, ymax = ci_upper),
      fill = "navy", 
      alpha = 0.2
    ) +
    # Spline line
    geom_line(
      data = express_partial,
      aes(x = distance, y = pred_flow, color = "Spline"),
      linewidth = 1.2
    ) +
    # Exponential line (dashed)
    geom_line(
      data = exp_results$predictions,
      aes(x = distance, y = pred_flow_exp, color = "Exponential"),
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    scale_color_manual(
      name = "Model",
      values = c("Spline" = "navy", "Exponential" = "darkred")
    ) +
    labs(
      x = "Distance (km)",
      y = "Flow"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6)
    )
} else {
  p1 <- ggplot(express_partial, aes(x = distance)) +
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper),
      fill = "navy", 
      alpha = 0.2
    ) +
    geom_line(
      aes(y = pred_flow),
      linewidth = 1.2,
      color = "navy"
    ) +
    labs(
      x = "Distance (km)",
      y = "Flow"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}

print(p1)

# --- Plot 2: Scaled version with BOTH spline and exponential curves ---
express_partial_scaled <- express_partial %>%
  mutate(
    flow_sum = sum(pred_flow),
    pred_flow_scaled = pred_flow / flow_sum,
    ci_lower_scaled = ci_lower / flow_sum,
    ci_upper_scaled = ci_upper / flow_sum,
    model_type = "Spline"
  )

# Scale exponential predictions
if (!is.null(exp_results)) {
  exp_partial_scaled <- exp_results$predictions %>%
    mutate(
      flow_sum = sum(pred_flow_exp),
      pred_flow_scaled = pred_flow_exp / flow_sum,
      model_type = "Exponential"
    )
  
  p2 <- ggplot() +
    # Spline confidence interval
    geom_ribbon(
      data = express_partial_scaled,
      aes(x = distance, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
      fill = "navy",
      alpha = 0.2
    ) +
    # Spline line
    geom_line(
      data = express_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Spline"),
      linewidth = 1.2
    ) +
    # Exponential line (dashed)
    geom_line(
      data = exp_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Exponential"),
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    scale_color_manual(
      name = "Model",
      values = c("Spline" = "navy", "Exponential" = "darkred")
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    ) +
    coord_cartesian(ylim = c(0.0013, 0.042), xlim = c(0.6, 12.1)) +
    labs(
      x = "Distance (km)",
      y = "Share of Flow (% of total)"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6)
    )
} else {
  # If exponential fit failed, just plot spline
  p2 <- ggplot(express_partial_scaled, aes(x = distance)) +
    geom_ribbon(
      aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
      fill = "navy",
      alpha = 0.2
    ) +
    geom_line(
      aes(y = pred_flow_scaled),
      linewidth = 1.2,
      color = "navy"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    ) +
    coord_cartesian(ylim = c(0.0013, 0.042), xlim = c(0.6, 12.1)) +
    labs(
      x = "Distance (km)",
      y = "Share of Flow (% of total)"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13)
    )
}

print(p2)

# --- Print detailed model summary ---
cat("\n\n=== DETAILED SPLINE MODEL SUMMARY ===\n")
cat(rep("=", 60), "\n", sep = "")
print(summary(results$model))

if (!is.null(exp_results)) {
  cat("\n\n=== DETAILED EXPONENTIAL MODEL SUMMARY ===\n")
  cat(rep("=", 60), "\n", sep = "")
  print(summary(exp_results$model))
}

# The exponential model was fit to unscaled data
exp_coefs <- coef(exp_results$model)

# For scaled data, only 'a' changes (gets divided by flow_sum)
flow_sum <- sum(express_partial$pred_flow)
a_scaled <- exp_coefs["a"] / flow_sum
b_scaled <- exp_coefs["b"]  # b stays the same!

cat("\n=== EXPONENTIAL PARAMETERS (SCALED) ===\n")
cat("a (scaled):", round(a_scaled, 6), "\n")
cat("b (decay rate):", round(b_scaled, 4), "\n")
