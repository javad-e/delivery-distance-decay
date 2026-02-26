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


# NEW FUNCTION: Fit exponential decay model to partial effect predictions
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


# NEW FUNCTION: Fit power law model to partial effect predictions
fit_power_law_model <- function(partial_effect_data, dist_min = 0, dist_max = 18, pred_points = 200) {
  # Fit power law to the spline partial effect predictions
  # flow = a * distance^(-b)
  
  # Filter out very small distances
  data_filtered <- partial_effect_data %>%
    filter(distance >= 0.5)  # Start from 0.5km to avoid near-zero issues
  
  cat("\n=== POWER LAW MODEL ===\n")
  cat("Fitting to", nrow(data_filtered), "data points\n")
  cat("Distance range:", min(data_filtered$distance), "to", max(data_filtered$distance), "km\n\n")
  
  # Method 1: Weighted log-linear regression (gives more weight to fit across all distances)
  data_filtered <- data_filtered %>%
    mutate(
      log_distance = log(distance),
      log_flow = log(pmax(pred_flow, 1))
    )
  
  # Fit weighted log-linear model
  # Use inverse variance weighting or uniform weights
  lm_fit <- lm(log_flow ~ log_distance, data = data_filtered)
  
  # Extract coefficients
  log_a <- coef(lm_fit)[1]
  b_coefficient <- coef(lm_fit)[2]
  
  # For power law: flow = a * distance^(-b)
  # In log space: log(flow) = log(a) - b*log(distance)
  # So b_coefficient from regression = -b, therefore b = -b_coefficient
  start_a <- exp(log_a)
  start_b <- -b_coefficient  # No abs() - keep the sign
  
  cat("Log-linear regression results:\n")
  cat("  Intercept (log_a):", log_a, "\n")
  cat("  Slope (b_coefficient):", b_coefficient, "\n")
  cat("  Implied a:", start_a, "\n")
  cat("  Implied b:", start_b, "\n\n")
  
  # If b is negative, the relationship is actually increasing, not decaying
  # This suggests power law may not be appropriate
  if (start_b < 0) {
    cat("WARNING: Implied b is negative, suggesting increasing relationship.\n")
    cat("Power law decay may not be appropriate for this data.\n\n")
  }
  
  tryCatch({
    # Method 2: Try NLS with multiple starting values
    best_model <- NULL
    best_aic <- Inf
    
    # Try different starting values
    start_values_list <- list(
      list(a = start_a, b = max(start_b, 0.1)),
      list(a = start_a * 1.5, b = max(start_b, 0.1) * 1.5),
      list(a = start_a * 0.5, b = max(start_b, 0.1) * 0.5),
      list(a = max(data_filtered$pred_flow) * 2, b = 1.0),
      list(a = max(data_filtered$pred_flow) * 2, b = 0.5)
    )
    
    for (i in seq_along(start_values_list)) {
      start_vals <- start_values_list[[i]]
      cat("Trying starting values: a =", start_vals$a, ", b =", start_vals$b, "\n")
      
      try_model <- try({
        nls(
          pred_flow ~ a * distance^(-b),
          data = data_filtered,
          start = start_vals,
          control = nls.control(maxiter = 1000, warnOnly = TRUE)
        )
      }, silent = TRUE)
      
      if (!inherits(try_model, "try-error")) {
        aic_val <- AIC(try_model)
        cat("  Success! AIC =", aic_val, "\n")
        if (aic_val < best_aic) {
          best_model <- try_model
          best_aic <- aic_val
        }
      } else {
        cat("  Failed\n")
      }
    }
    
    if (!is.null(best_model)) {
      model_power <- best_model
      
      # Generate predictions
      distance_grid <- tibble(distance = seq(0.5, dist_max, length.out = pred_points))
      distance_grid$pred_flow_power <- predict(model_power, newdata = distance_grid)
      
      # Calculate R² (fit to partial effect curve)
      # Use only the filtered distance range for fair comparison
      comparison_data <- partial_effect_data %>%
        filter(distance >= 0.5)
      
      predictions_at_data <- predict(model_power, newdata = comparison_data)
      ss_total <- sum((comparison_data$pred_flow - mean(comparison_data$pred_flow))^2)
      ss_res <- sum((comparison_data$pred_flow - predictions_at_data)^2)
      r2_power <- 1 - ss_res / ss_total
      
      cat("\n=== BEST FIT RESULTS ===\n")
      cat("Formula: flow = a * distance^(-b)\n")
      cat("Parameters:\n")
      print(coef(model_power))
      cat("R²:", r2_power, "\n")
      cat("AIC:", best_aic, "\n\n")
      
      list(
        predictions = distance_grid %>%
          mutate(r2 = r2_power),
        model = model_power,
        r2 = r2_power
      )
    } else {
      stop("All NLS attempts failed")
    }
    
  }, error = function(e) {
    cat("\nAll NLS fitting attempts failed, using log-linear approximation\n")
    cat("Error was:", e$message, "\n\n")
    
    # Fallback: use the log-linear fit
    distance_grid <- tibble(distance = seq(0.5, dist_max, length.out = pred_points))
    distance_grid$log_distance <- log(distance_grid$distance)
    distance_grid$pred_flow_power <- exp(predict(lm_fit, newdata = distance_grid))
    
    # Calculate R² on original scale
    comparison_data <- partial_effect_data %>%
      filter(distance >= 0.5)
    comparison_data$log_distance <- log(comparison_data$distance)
    
    pred_fitted <- exp(predict(lm_fit, newdata = comparison_data))
    ss_total <- sum((comparison_data$pred_flow - mean(comparison_data$pred_flow))^2)
    ss_res <- sum((comparison_data$pred_flow - pred_fitted)^2)
    r2_power <- 1 - ss_res / ss_total
    
    cat("Using log-linear approximation\n")
    cat("Formula: flow = exp(", log_a, ") * distance^(", b_coefficient, ")\n")
    cat("Equivalent: flow =", start_a, "* distance^", b_coefficient, "\n")
    cat("R² (on original scale):", r2_power, "\n\n")
    
    list(
      predictions = distance_grid %>%
        select(distance, pred_flow_power) %>%
        mutate(r2 = r2_power),
      model = lm_fit,
      r2 = r2_power
    )
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

# Fit power law model to the partial effect curve
power_results <- fit_power_law_model(express_partial, dist_max = 18)

# --- Display Model Summary ---
cat("\n\n=== MODEL SUMMARY ===\n")
cat("Control variables used:", paste(results$control_vars_used, collapse = ", "), "\n")
cat("Spline Model R² (Conditional):", unique(express_partial$r2), "\n")
if (!is.null(exp_results)) {
  cat("Exponential Model R²:", exp_results$r2, "\n")
}
if (!is.null(power_results)) {
  cat("Power Law Model R²:", power_results$r2, "\n")
}

# --- Create 4 separate plots in a grid ---

# Plot 1: Unscaled with Power Law (first row, left)
if (!is.null(power_results)) {
  p1 <- ggplot() +
    geom_ribbon(
      data = express_partial,
      aes(x = distance, ymin = ci_lower, ymax = ci_upper),
      fill = "navy", 
      alpha = 0.2
    ) +
    geom_line(
      data = express_partial,
      aes(x = distance, y = pred_flow, color = "Spline"),
      linewidth = 1.2
    ) +
    geom_line(
      data = power_results$predictions,
      aes(x = distance, y = pred_flow_power, color = "Power Law"),
      linewidth = 1,
      linetype = "dotted"
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Spline" = "navy", "Power Law" = "#2E8B57")
    ) +
    scale_y_continuous(
      breaks = seq(0, 750, by = 100)
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    ) +
    coord_cartesian(xlim = c(0.6, 12.1), ylim = c(20, 700)) +
    labs(
      title = "(I) Power Law - Unscaled",
      x = "Distance (km)",
      y = "Flow"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6),
      legend.text = element_text(size = 11)
    )
} else {
  p1 <- ggplot(express_partial, aes(x = distance)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "navy", alpha = 0.2) +
    geom_line(aes(y = pred_flow), linewidth = 1.2, color = "navy") +
    scale_y_continuous(breaks = seq(0, 750, by = 100)) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
    coord_cartesian(xlim = c(0.6, 12.1), ylim = c(20, 700)) +
    labs(title = "(I) Power Law - Unscaled", x = "Distance (km)", y = "Flow") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13)
    )
}

# Plot 2: Scaled with Power Law (first row, right)
# Prepare scaled data
express_partial_scaled <- express_partial %>%
  mutate(
    flow_sum = sum(pred_flow),
    pred_flow_scaled = pred_flow / flow_sum,
    ci_lower_scaled = ci_lower / flow_sum,
    ci_upper_scaled = ci_upper / flow_sum
  )

if (!is.null(power_results)) {
  power_partial_scaled <- power_results$predictions %>%
    mutate(
      flow_sum = sum(pred_flow_power),
      pred_flow_scaled = pred_flow_power / flow_sum
    )
  
  p2 <- ggplot() +
    geom_ribbon(
      data = express_partial_scaled,
      aes(x = distance, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
      fill = "navy",
      alpha = 0.2
    ) +
    geom_line(
      data = express_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Spline"),
      linewidth = 1.2
    ) +
    geom_line(
      data = power_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Power Law"),
      linewidth = 1,
      linetype = "dotted"
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Spline" = "navy", "Power Law" = "#2E8B57")
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
      title = "(II) Power Law - Scaled",
      x = "Distance (km)",
      y = "Share of Flow (% of total)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6),
      legend.text = element_text(size = 11)
    )
} else {
  p2 <- ggplot(express_partial_scaled, aes(x = distance)) +
    geom_ribbon(aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled), fill = "navy", alpha = 0.2) +
    geom_line(aes(y = pred_flow_scaled), linewidth = 1.2, color = "navy") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Distance (km)", y = "Share of Flow (% of total)") +
    theme_minimal() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13))
}

# Plot 3: Unscaled with Exponential (second row, left)
if (!is.null(exp_results)) {
  p3 <- ggplot() +
    geom_ribbon(
      data = express_partial,
      aes(x = distance, ymin = ci_lower, ymax = ci_upper),
      fill = "navy", 
      alpha = 0.2
    ) +
    geom_line(
      data = express_partial,
      aes(x = distance, y = pred_flow, color = "Spline"),
      linewidth = 1.2
    ) +
    geom_line(
      data = exp_results$predictions,
      aes(x = distance, y = pred_flow_exp, color = "Exponential"),
      linewidth = 1,
      linetype = "longdash"
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Spline" = "navy", "Exponential" = "#B22222")
    ) +
    scale_y_continuous(
      breaks = seq(0, 750, by = 100)
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    ) +
    coord_cartesian(xlim = c(0.6, 12.1), ylim = c(20, 700)) +
    labs(
      title = "(III) Exponential - Unscaled",
      x = "Distance (km)",
      y = "Flow"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6),
      legend.text = element_text(size = 11)
    )
} else {
  p3 <- ggplot(express_partial, aes(x = distance)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "navy", alpha = 0.2) +
    geom_line(aes(y = pred_flow), linewidth = 1.2, color = "navy") +
    scale_y_continuous(breaks = seq(0, 750, by = 100)) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
    coord_cartesian(xlim = c(0.6, 12.1), ylim = c(20, 700)) +
    labs(title = "(III) Exponential - Unscaled", x = "Distance (km)", y = "Flow") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13)
    )
}

# Plot 4: Scaled with Exponential (second row, right)
if (!is.null(exp_results)) {
  exp_partial_scaled <- exp_results$predictions %>%
    mutate(
      flow_sum = sum(pred_flow_exp),
      pred_flow_scaled = pred_flow_exp / flow_sum
    )
  
  p4 <- ggplot() +
    geom_ribbon(
      data = express_partial_scaled,
      aes(x = distance, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
      fill = "navy",
      alpha = 0.2
    ) +
    geom_line(
      data = express_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Spline"),
      linewidth = 1.2
    ) +
    geom_line(
      data = exp_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Exponential"),
      linewidth = 1,
      linetype = "longdash"
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Spline" = "navy", "Exponential" = "#B22222")
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
      title = "(IV) Exponential - Scaled",
      x = "Distance (km)",
      y = "Share of Flow (% of total)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      legend.position = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin = margin(6, 6, 6, 6),
      legend.text = element_text(size = 11)
    )
} else {
  p4 <- ggplot(express_partial_scaled, aes(x = distance)) +
    geom_ribbon(aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled), fill = "navy", alpha = 0.2) +
    geom_line(aes(y = pred_flow_scaled), linewidth = 1.2, color = "navy") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Distance (km)", y = "Share of Flow (% of total)") +
    theme_minimal() +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13))
}

# Combine plots in a 2x2 grid: power law in first row, exponential in second row
library(patchwork)

combined_plot <- (p1 | p2) / (p3 | p4)

print(combined_plot)

# --- Print detailed model summary ---
cat("\n\n=== DETAILED SPLINE MODEL SUMMARY ===\n")
cat(rep("=", 60), "\n", sep = "")
print(summary(results$model))

if (!is.null(exp_results)) {
  cat("\n\n=== DETAILED EXPONENTIAL MODEL SUMMARY ===\n")
  cat(rep("=", 60), "\n", sep = "")
  print(summary(exp_results$model))
}

if (!is.null(power_results)) {
  cat("\n\n=== DETAILED POWER LAW MODEL SUMMARY ===\n")
  cat(rep("=", 60), "\n", sep = "")
  print(summary(power_results$model))
}



# --- New Figure: Scaled Spline with Both Power Law and Exponential Overlaid ---
# --- New Figure: Scaled Spline with Both Power Law and Exponential Overlaid ---

if (!is.null(exp_results) && !is.null(power_results)) {
  
  exp_partial_scaled <- exp_results$predictions %>%
    mutate(
      flow_sum = sum(pred_flow_exp),
      pred_flow_scaled = pred_flow_exp / flow_sum
    )
  
  power_partial_scaled <- power_results$predictions %>%
    mutate(
      flow_sum = sum(pred_flow_power),
      pred_flow_scaled = pred_flow_power / flow_sum
    )
  
  p_combined <- ggplot() +
    # Confidence ribbon for spline
    geom_ribbon(
      data = express_partial_scaled,
      aes(x = distance, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
      fill = "navy",
      alpha = 0.15
    ) +
    # Spline curve
    geom_line(
      data = express_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Spline", linetype = "Spline"),
      linewidth = 1.3
    ) +
    # Power Law curve
    geom_line(
      data = power_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Power Law", linetype = "Power Law"),
      linewidth = 0.95
    ) +
    # Exponential curve
    geom_line(
      data = exp_partial_scaled,
      aes(x = distance, y = pred_flow_scaled, color = "Exponential", linetype = "Exponential"),
      linewidth = 0.95
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Spline"      = "navy",
        "Power Law"   = "#2E8B57",   # firebrick — deep red, print-safe
        "Exponential" = "#B22222"    # sea green — muted, distinct from red 
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Spline"      = "solid",
        "Power Law"   = "dotted",  # clean long dash, readable at small sizes
        "Exponential" = "longdash"     # standard dash, clearly distinct
      )
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)
    ) +
    scale_x_continuous(
      breaks = 1:12
    ) +
    coord_cartesian(
      xlim = c(0.6, 12.1),
      ylim = c(0.0013, 0.042)
    ) +
    labs(
      x = "Distance (km)",
      y = "Share of Flow (% of total)"
    ) +
    theme_minimal() +
    theme(
      panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 13),
      legend.position      = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.margin        = margin(6, 6, 6, 6),
      legend.text          = element_text(size = 12),
      legend.key.width     = unit(1.8, "cm")   # widen key so dashes render clearly
    )
  
  print(p_combined)
  
} else {
  warning("One or both parametric models failed to fit; combined overlay plot skipped.")
}