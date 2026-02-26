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
  
  beta <- fixef(model)
  vcov_beta <- as.matrix(vcov(model))
  
  # Calculate predictions
  distance_grid$pred_log_flow <- as.vector(X_grid %*% beta)
  
  distance_grid$se_log_flow <- sqrt(diag(X_grid %*% vcov_beta %*% t(X_grid)))
  
  distance_grid <- distance_grid %>%
    mutate(
      log_ci_lower = pred_log_flow - 1.96 * se_log_flow,
      log_ci_upper = pred_log_flow + 1.96 * se_log_flow,
      pred_flow = pmax(exp(pred_log_flow) - 1, 0),
      ci_lower = pmax(exp(log_ci_lower) - 1, 0),
      ci_upper = pmax(exp(log_ci_upper) - 1, 0)
    )
  
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

# --- Process Dubai Data ---
cat("\n=== Processing Dubai Dataset ===\n")
df <- read_csv("~/imperial/Reach of Last Mile/platform1_dubai_curve_estimation.csv", show_col_types = FALSE)
cat("Dataset dimensions:", nrow(df), "rows,", ncol(df), "columns\n")

results <- run_spline_model_return_curve(df, 
                                         control_vars = CONTROL_VARIABLES,
                                         dist_max = 18, 
                                         spline_step = 2)

express_partial <- results$partial_effect
express_raw <- results$raw_data

# --- Display Model Summary ---
cat("\n\n=== MODEL SUMMARY ===\n")
cat("Control variables used:", paste(results$control_vars_used, collapse = ", "), "\n")
cat("Model R² (Conditional):", unique(express_partial$r2), "\n")

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
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(p1)

express_partial_scaled <- express_partial %>%
  mutate(
    # Calculate the sum of all predicted flows
    flow_sum = sum(pred_flow),
    # Scale prediction AND CI bounds by dividing by the sum
    pred_flow_scaled = pred_flow / flow_sum,
    ci_lower_scaled = ci_lower / flow_sum,
    ci_upper_scaled = ci_upper / flow_sum
  )

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
    breaks=c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04)

  ) +
  scale_x_continuous(
    breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),
  ) +
  coord_cartesian(ylim = c(0.0013, 0.042), xlim = c(0.6,12.1)) +
  labs(
    x = "Distance (km)",
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

print(p2)

saveRDS(express_partial_scaled, "~/imperial/Reach of Last Mile/platform1_curves.rds")

# --- Print detailed model summary ---
cat("\n\n=== DETAILED MODEL SUMMARY ===\n")
cat(rep("=", 60), "\n", sep = "")
print(summary(results$model))
