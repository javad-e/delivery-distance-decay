library(tidyverse)
library(mgcv)
library(scales)

# Define control variables to experiment with
CONTROL_VARIABLES <- c(
  "num_vendors_dest",
  "income_cat_origin",
  "income_cat_dest",
  "vendor_entropy_dest"
)


run_gp_model_return_curve <- function(df, 
                                      control_vars = CONTROL_VARIABLES,
                                      dist_min = 0, 
                                      dist_max = 18, 
                                      pred_points = 160,
                                      n_basis = 30,
                                      include_random_effects = TRUE) {  # NEW PARAMETER
  
  # Filter and prepare data
  df <- df %>%
    filter(euclidean_km >= dist_min, euclidean_km <= dist_max) %>%
    mutate(
      log_flow = log(flow + 1),
      origin = factor(origin),           # ADD THIS
      destination = factor(destination)  # ADD THIS
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
  
  # Build formula with GP term and linear control variables
  control_terms <- ""
  if (length(available_controls) > 0) {
    control_terms_std <- paste0(available_controls, "_std")
    control_terms <- paste0(" + ", paste(control_terms_std, collapse = " + "))
  }
  
  # Add random effects if requested
  random_effects <- ""
  if (include_random_effects) {
    random_effects <- " + s(origin, bs='re') + s(destination, bs='re')"
  }
  
  # Using s() with bs="gp" for Gaussian Process
  formula_text <- paste0("log_flow ~ 1 + s(euclidean_km, bs='gp', k=", n_basis, 
                         ", m=c(3, 0.5))", control_terms, random_effects)
  
  cat("Model formula:", formula_text, "\n")
  cat("GP kernel: Matérn 5/2 with", n_basis, "basis functions\n")
  cat("Random effects:", ifelse(include_random_effects, "YES (origin + destination)", "NO"), "\n\n")
  
  # Fit GAM with GP smooth
  model <- gam(as.formula(formula_text), data = df, method = "REML")
  
  # Generate predictions on grid
  distance_grid <- tibble(euclidean_km = seq(dist_min, dist_max, length.out = pred_points))
  
  # Set control variables to their mean (0 for standardized variables)
  for (var in available_controls) {
    distance_grid[[paste0(var, "_std")]] <- 0
  }
  
  # Set random effects to a reference level (for marginal prediction)
  if (include_random_effects) {
    distance_grid$origin <- factor(levels(df$origin)[1], levels = levels(df$origin))
    distance_grid$destination <- factor(levels(df$destination)[1], levels = levels(df$destination))
  }
  
  # Get predictions with standard errors
  # Use exclude to exclude random effects from predictions (like your spline model does)
  if (include_random_effects) {
    preds <- predict(model, newdata = distance_grid, se.fit = TRUE, 
                     exclude = c("s(origin)", "s(destination)"),
                     unconditional = FALSE)
  } else {
    preds <- predict(model, newdata = distance_grid, se.fit = TRUE, 
                     unconditional = FALSE)
  }
  
  # Calculate confidence intervals on log scale, then transform
  distance_grid <- distance_grid %>%
    mutate(
      pred_log_flow = preds$fit,
      se_log_flow = preds$se.fit,
      log_ci_lower = pred_log_flow - 1.8 * se_log_flow,
      log_ci_upper = pred_log_flow + 1.8 * se_log_flow,
      pred_flow = pmax(exp(pred_log_flow) - 1, 0),
      ci_lower = pmax(exp(log_ci_lower) - 1, 0),
      ci_upper = pmax(exp(log_ci_upper) - 1, 0)
    )
  
  # Calculate R²
  df$pred_log_flow <- predict(model)
  ss_total <- sum((df$log_flow - mean(df$log_flow))^2)
  ss_res <- sum((df$log_flow - df$pred_log_flow)^2)
  r2 <- 1 - ss_res / ss_total
  
  # Extract GP-specific information
  smooth_info <- model$smooth[[1]]
  length_scale <- sqrt(1 / smooth_info$delta)
  
  # Calculate derivative of GP function (rate of decay)
  eps <- 0.01
  distance_grid_plus <- distance_grid %>% mutate(euclidean_km = euclidean_km + eps)
  distance_grid_minus <- distance_grid %>% mutate(euclidean_km = euclidean_km - eps)
  
  if (include_random_effects) {
    pred_plus <- predict(model, newdata = distance_grid_plus, se.fit = TRUE,
                         exclude = c("s(origin)", "s(destination)"),
                         unconditional = FALSE)
    pred_minus <- predict(model, newdata = distance_grid_minus, se.fit = TRUE,
                          exclude = c("s(origin)", "s(destination)"),
                          unconditional = FALSE)
  } else {
    pred_plus <- predict(model, newdata = distance_grid_plus, se.fit = TRUE,
                         unconditional = FALSE)
    pred_minus <- predict(model, newdata = distance_grid_minus, se.fit = TRUE,
                          unconditional = FALSE)
  }
  
  distance_grid <- distance_grid %>%
    mutate(
      derivative = (pred_plus$fit - pred_minus$fit) / (2 * eps),
      derivative_se = sqrt(pred_plus$se.fit^2 + pred_minus$se.fit^2) / (2 * eps),
      derivative_ci_lower = derivative - 1.96 * derivative_se,
      derivative_ci_upper = derivative + 1.96 * derivative_se
    )
  
  # Return results
  list(
    partial_effect = distance_grid %>%
      dplyr::select(distance = euclidean_km, pred_flow, ci_lower, ci_upper,
                    pred_log_flow, se_log_flow,
                    derivative, derivative_se, derivative_ci_lower, derivative_ci_upper) %>%
      mutate(r2 = r2),
    raw_data = df %>%
      dplyr::select(distance = euclidean_km, flow),
    model = model,
    control_vars_used = available_controls,
    gp_info = list(
      length_scale = length_scale,
      edf = sum(model$edf),
      n_basis = n_basis,
      has_random_effects = include_random_effects
    )
  )
}

# --- Process Dubai Data ---
cat("\n=== Processing Dubai Dataset ===\n")
df <- read_csv("~/imperial/Reach of Last Mile/platform1_dubai_curve_estimation.csv", 
               show_col_types = FALSE)
cat("Dataset dimensions:", nrow(df), "rows,", ncol(df), "columns\n")

# NOW WITH RANDOM EFFECTS (comparable to spline model)
results <- run_gp_model_return_curve(df, 
                                     control_vars = CONTROL_VARIABLES,
                                     dist_max = 18,
                                     n_basis = 30,
                                     include_random_effects = TRUE)  # KEY CHANGE

# Extract results
express_partial <- results$partial_effect
express_raw <- results$raw_data

# --- Display Model Summary ---
cat("\n\n=== MODEL SUMMARY ===\n")
cat("Control variables used:", paste(results$control_vars_used, collapse = ", "), "\n")
cat("Model R²:", unique(express_partial$r2), "\n")
cat("GP Length Scale:", round(results$gp_info$length_scale, 3), "km\n")
cat("Effective degrees of freedom:", round(results$gp_info$edf, 2), "\n")
cat("Number of basis functions:", results$gp_info$n_basis, "\n")
cat("Random effects:", ifelse(results$gp_info$has_random_effects, "YES", "NO"), "\n")

# --- Plot 1: Partial Effect Curve with 95% CI (Original Scale) ---
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
    y = "Predicted Flow"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

print(p1)

# --- Plot 2: Scaled to [0,1] ---
express_partial_scaled <- express_partial %>%
  mutate(
    flow_sum = sum(pred_flow),
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
  coord_cartesian(ylim = c(0.0013, 0.042), xlim = c(0.6,12)) +
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