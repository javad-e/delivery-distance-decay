library(tidyverse)
library(lme4)
library(splines)
library(scales)

# ===========================
# CONFIGURABLE PARAMETERS
# ===========================

# Upper bounds for each metric
EUCLIDEAN_MAX <- 18
H3_EUCLIDEAN_MAX <- 18
DRIVING_MAX <- 45
TRAVEL_TIME_MAX <- 90  # in minutes

# Spline step size for each metric
EUCLIDEAN_STEP <- 2
H3_EUCLIDEAN_STEP <- 2
DRIVING_STEP <- 4
TRAVEL_TIME_STEP <- 10

# Plot parameters for scaled plots - Euclidean
EUCLIDEAN_YLIM <- c(0.0013, 0.035)
EUCLIDEAN_XLIM <- c(0.6, 12.5)
EUCLIDEAN_X_BREAKS <- c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5)
EUCLIDEAN_Y_BREAKS <- c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)

# Plot parameters for scaled plots - H3 Euclidean
H3_EUCLIDEAN_YLIM <- c(0.0013, 0.035)
H3_EUCLIDEAN_XLIM <- c(0.6, 12.5)
H3_EUCLIDEAN_X_BREAKS <- c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5)
H3_EUCLIDEAN_Y_BREAKS <- c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)

# Plot parameters for scaled plots - Driving
DRIVING_YLIM <- c(0.0013, 0.035)
DRIVING_XLIM <- c(1.3, 28)
DRIVING_X_BREAKS <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
DRIVING_Y_BREAKS <- c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)

# Plot parameters for scaled plots - Travel Time
TRAVEL_TIME_YLIM <- c(0.0013, 0.035)
TRAVEL_TIME_XLIM <- c(3, 55)
TRAVEL_TIME_X_BREAKS <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55)
TRAVEL_TIME_Y_BREAKS <- c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)

# Add or remove column names from this list to change the model specification
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

# ===========================
# MAIN FUNCTION
# ===========================

run_spline_model_return_curve <- function(df, 
                                          distance_col,
                                          control_vars = CONTROL_VARIABLES,
                                          dist_min = 0, 
                                          dist_max = 15, 
                                          spline_step = 3, 
                                          pred_points = 200) {
  
  # Check if distance column exists
  if (!distance_col %in% names(df)) {
    stop(paste("Distance column", distance_col, "not found in dataset"))
  }
  
  # Filter and prepare data
  df <- df %>%
    filter(.data[[distance_col]] >= dist_min, .data[[distance_col]] <= dist_max) %>%
    mutate(
      origin = factor(origin),
      destination = factor(destination),
      log_flow = log(flow + 1),
      distance_var = .data[[distance_col]]  # Create a generic distance variable
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
    df$distance_var,
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
  
  # Set control variables to their mean (0 for standardized variables)
  for (var in available_controls) {
    distance_grid[[paste0(var, "_std")]] <- 0
  }
  
  # Get design matrix for fixed effects predictions
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
      dplyr::select(distance = distance_var, flow),
    model = model,
    control_vars_used = available_controls
  )
}

# ===========================
# LOAD DATA
# ===========================
cat("\n=== Loading Dubai Dataset ===\n")
df <- read_csv("~/imperial/Reach of Last Mile/platform1_dubai_curve_estimation.csv", show_col_types = FALSE)
cat("Dataset dimensions:", nrow(df), "rows,", ncol(df), "columns\n")

# ===========================
# RUN ANALYSES FOR EACH METRIC
# ===========================

metrics <- list(
  list(
    col = "euclidean_km", 
    max = EUCLIDEAN_MAX, 
    step = EUCLIDEAN_STEP, 
    label = "Euclidean Distance",
    ylim = EUCLIDEAN_YLIM,
    xlim = EUCLIDEAN_XLIM,
    x_breaks = EUCLIDEAN_X_BREAKS,
    y_breaks = EUCLIDEAN_Y_BREAKS
  ),
  list(
    col = "h3_euclidean_km", 
    max = H3_EUCLIDEAN_MAX, 
    step = H3_EUCLIDEAN_STEP, 
    label = "H3 Euclidean Distance",
    ylim = H3_EUCLIDEAN_YLIM,
    xlim = H3_EUCLIDEAN_XLIM,
    x_breaks = H3_EUCLIDEAN_X_BREAKS,
    y_breaks = H3_EUCLIDEAN_Y_BREAKS
  ),
  list(
    col = "driving_km", 
    max = DRIVING_MAX, 
    step = DRIVING_STEP, 
    label = "Driving Distance",
    ylim = DRIVING_YLIM,
    xlim = DRIVING_XLIM,
    x_breaks = DRIVING_X_BREAKS,
    y_breaks = DRIVING_Y_BREAKS
  ),
  list(
    col = "travel_time", 
    max = TRAVEL_TIME_MAX, 
    step = TRAVEL_TIME_STEP, 
    label = "Travel Time",
    ylim = TRAVEL_TIME_YLIM,
    xlim = TRAVEL_TIME_XLIM,
    x_breaks = TRAVEL_TIME_X_BREAKS,
    y_breaks = TRAVEL_TIME_Y_BREAKS
  )
)

results_list <- list()
plots_original <- list()
plots_scaled <- list()

for (metric in metrics) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("=== Processing", metric$label, "(", metric$col, ") ===\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Check if column exists
  if (!metric$col %in% names(df)) {
    cat("WARNING: Column", metric$col, "not found in dataset. Skipping...\n")
    next
  }
  
  # Print basic statistics about this metric
  cat("\nData diagnostics for", metric$col, ":\n")
  cat("  Min:", min(df[[metric$col]], na.rm = TRUE), "\n")
  cat("  Max:", max(df[[metric$col]], na.rm = TRUE), "\n")
  cat("  Mean:", mean(df[[metric$col]], na.rm = TRUE), "\n")
  cat("  Median:", median(df[[metric$col]], na.rm = TRUE), "\n")
  cat("  NA values:", sum(is.na(df[[metric$col]])), "\n")
  cat("  Unique values:", length(unique(df[[metric$col]])), "\n")
  
  # Show distribution
  cat("\nQuantiles:\n")
  print(quantile(df[[metric$col]], probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE))
  
  # Run model
  results <- run_spline_model_return_curve(
    df, 
    distance_col = metric$col,
    control_vars = CONTROL_VARIABLES,
    dist_max = metric$max, 
    spline_step = metric$step
  )
  
  # Store results
  results_list[[metric$col]] <- results
  
  # Extract results
  partial_effect <- results$partial_effect
  raw_data <- results$raw_data
  
  # Display model summary
  cat("\n=== MODEL SUMMARY ===\n")
  cat("Control variables used:", paste(results$control_vars_used, collapse = ", "), "\n")
  cat("Model R² (Conditional):", unique(partial_effect$r2), "\n\n")
  
  # ===========================
  # PLOT 1: Original Scale
  # ===========================
  
  x_label <- if (metric$col == "travel_time") "Travel Time (minutes)" else "Distance (km)"
  
  p1 <- ggplot(partial_effect, aes(x = distance)) +
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
      title = paste("Partial Effect of", metric$label, "on Flow (Dubai)"),
      x = x_label,
      y = "Predicted Flow"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  plots_original[[metric$col]] <- p1
  print(p1)
  
  # ===========================
  # PLOT 2: Scaled (% of total)
  # ===========================
  
  partial_effect_scaled <- partial_effect %>%
    mutate(
      flow_sum = sum(pred_flow),
      pred_flow_scaled = pred_flow / flow_sum,
      ci_lower_scaled = ci_lower / flow_sum,
      ci_upper_scaled = ci_upper / flow_sum
    )
  
  p2 <- ggplot(partial_effect_scaled, aes(x = distance)) +
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
      breaks = metric$y_breaks
    ) +
    scale_x_continuous(
      breaks = metric$x_breaks
    ) +
    coord_cartesian(
      ylim = metric$ylim, 
      xlim = metric$xlim
    ) +
    labs(
      title = paste(metric$label),
      x = x_label,
      y = "Share of Flow (% of total)"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13)
    )
  
  plots_scaled[[metric$col]] <- p2
  print(p2)
  
  # Save individual results
  saveRDS(partial_effect_scaled, 
          paste0("~/imperial/Reach of Last Mile/platform1_curve_", metric$col, ".rds"))
  
  cat("\n")
}

# ===========================
# SAVE ALL RESULTS
# ===========================

saveRDS(results_list, "~/imperial/Reach of Last Mile/platform1_curves_all_metrics.rds")
cat("\n=== All results saved ===\n")
cat("Individual curves saved as: platform1_curve_[metric].rds\n")
cat("Combined results saved as: platform1_curves_all_metrics.rds\n")

# ===========================
# PRINT DETAILED MODEL SUMMARIES
# ===========================

for (metric in metrics) {
  if (metric$col %in% names(results_list)) {
    cat("\n\n", rep("=", 70), "\n", sep = "")
    cat("=== DETAILED MODEL SUMMARY:", metric$label, "===\n")
    cat(rep("=", 70), "\n", sep = "")
    print(summary(results_list[[metric$col]]$model))
  }
}