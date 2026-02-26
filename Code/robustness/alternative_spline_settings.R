library(tidyverse)
library(lme4)
library(mgcv)  # Added for penalized splines
library(splines)
library(scales)

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


run_spline_model_return_curve <- function(df, 
                                          control_vars = CONTROL_VARIABLES,
                                          dist_min = 0, 
                                          dist_max = 18, 
                                          spline_step = 3, 
                                          pred_points = 200,
                                          spline_type = "bs",
                                          spline_df = NULL,
                                          spline_degree = 3,
                                          spline_k = 10) {  # Added k parameter
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
  
  # Handle penalized splines separately with GAM
  if (spline_type == "ps") {
    # Build formula for GAM with penalized spline (cubic P-spline by default)
    control_terms <- ""
    if (length(available_controls) > 0) {
      control_terms_std <- paste0(available_controls, "_std")
      control_terms <- paste0(" + ", paste(control_terms_std, collapse = " + "))
    }
    
    # bs='ps' uses cubic P-splines by default (m=c(2,1) which is cubic with first-order penalty)
    formula_text <- paste0("log_flow ~ s(euclidean_km, bs='ps', k=", spline_k, ", m=c(2,1))", 
                           control_terms,
                           " + s(origin, bs='re') + s(destination, bs='re')")
    
    cat("Model formula:", formula_text, "\n\n")
    
    # Fit GAM model
    model <- gam(as.formula(formula_text), data = df, method = "REML")
    
    # Generate predictions on grid
    distance_grid <- tibble(distance = seq(dist_min, dist_max, length.out = pred_points))
    distance_grid$euclidean_km <- distance_grid$distance
    distance_grid$origin <- factor(levels(df$origin)[1], levels = levels(df$origin))
    distance_grid$destination <- factor(levels(df$destination)[1], levels = levels(df$destination))
    
    # Set control variables to their mean (0 for standardized variables)
    for (var in available_controls) {
      distance_grid[[paste0(var, "_std")]] <- 0
    }
    
    # Get predictions with standard errors
    preds <- predict(model, newdata = distance_grid, se.fit = TRUE, 
                     exclude = c("s(origin)", "s(destination)"))
    
    distance_grid$pred_log_flow <- preds$fit
    distance_grid$se_log_flow <- preds$se.fit
    
    # Calculate confidence intervals on log scale, then transform
    distance_grid <- distance_grid %>%
      mutate(
        log_ci_lower = pred_log_flow - 1.96 * se_log_flow,
        log_ci_upper = pred_log_flow + 1.96 * se_log_flow,
        pred_flow = pmax(exp(pred_log_flow) - 1, 0),
        ci_lower = pmax(exp(log_ci_lower) - 1, 0),
        ci_upper = pmax(exp(log_ci_upper) - 1, 0)
      )
    
    # Calculate R² equivalent for GAM
    r2_conditional <- summary(model)$r.sq
    model_aic <- AIC(model)
    model_bic <- BIC(model)
    
  } else {
    # Original lmer approach for bs, ns, poly
    # Create spline basis based on type
    if (spline_type == "bs") {
      # B-spline (cubic by default)
      knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
      spline_basis <- bs(
        df$euclidean_km,
        knots = knots,
        degree = spline_degree,
        intercept = FALSE,
        Boundary.knots = c(dist_min, dist_max)
      )
    } else if (spline_type == "ns") {
      # Natural cubic spline
      if (!is.null(spline_df)) {
        spline_basis <- ns(df$euclidean_km, df = spline_df)
      } else {
        knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
        spline_basis <- ns(df$euclidean_km, knots = knots)
      }
    } else if (spline_type == "poly") {
      # Polynomial
      spline_basis <- poly(df$euclidean_km, degree = spline_degree, raw = FALSE)
    }
    
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
    
    # Create matching spline basis for prediction
    if (spline_type == "bs") {
      spline_basis_grid <- bs(
        distance_grid$distance,
        knots = knots,
        degree = spline_degree,
        intercept = FALSE,
        Boundary.knots = c(dist_min, dist_max)
      )
    } else if (spline_type == "ns") {
      if (!is.null(spline_df)) {
        spline_basis_grid <- ns(distance_grid$distance, df = spline_df)
      } else {
        spline_basis_grid <- ns(distance_grid$distance, knots = knots)
      }
    } else if (spline_type == "poly") {
      spline_basis_grid <- predict(poly(df$euclidean_km, degree = spline_degree, raw = FALSE), 
                                   distance_grid$distance)
    }
    
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
    
    # Calculate AIC and BIC for model comparison
    model_aic <- AIC(model)
    model_bic <- BIC(model)
  }
  
  # Return partial effect curve and raw data with parameter uncertainty CIs
  list(
    partial_effect = distance_grid %>%
      dplyr::select(distance, pred_flow, ci_lower, ci_upper) %>%
      mutate(r2 = r2_conditional, aic = model_aic, bic = model_bic),
    raw_data = df %>%
      dplyr::select(distance = euclidean_km, flow),
    model = model,
    control_vars_used = available_controls,
    spline_type = spline_type,
    spline_params = list(step = spline_step, df = spline_df, degree = spline_degree, k = spline_k)
  )
}

# --- Process Dubai Data ---
cat("\n=== Processing Dubai Dataset ===\n")
df <- read_csv("~/imperial/Reach of Last Mile/platform1_dubai_curve_estimation.csv", show_col_types = FALSE)
cat("Dataset dimensions:", nrow(df), "rows,", ncol(df), "columns\n")

# Define robustness check specifications
# First row: Compare cubic natural splines, penalized splines (cubic), and B-splines (cubic)
# Subsequent rows: Other specifications
spline_specs <- list(
  # FIRST ROW: Spline type comparison (all cubic)
  list(name = "Natural Spline (cubic)", 
       type = "ns", step = 3, df = NULL, degree = NULL, k = 10),
  
  list(name = "Penalized Spline (cubic)", 
       type = "ps", step = NULL, df = NULL, degree = NULL, k = 10),
  
  list(name = "B-spline (cubic)", 
       type = "bs", step = 3, df = NULL, degree = 3, k = 10),
  
  # SUBSEQUENT ROWS: Other variations
  list(name = "Polynomial (degree 3)", 
       type = "poly", step = NULL, df = NULL, degree = 3, k = NULL),
  
  list(name = "Polynomial (degree 4)", 
       type = "poly", step = NULL, df = NULL, degree = 4, k = NULL),
  
  list(name = "Polynomial (degree 5)", 
       type = "poly", step = NULL, df = NULL, degree = 5, k = NULL)
)

# Run models for all specifications
all_results <- list()
model_comparison <- tibble()

for (i in seq_along(spline_specs)) {
  spec <- spline_specs[[i]]
  cat("\n\n", rep("=", 70), "\n", sep = "")
  cat("Running:", spec$name, "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  results <- run_spline_model_return_curve(
    df, 
    control_vars = CONTROL_VARIABLES,
    dist_max = 18, 
    spline_step = ifelse(is.null(spec$step), 3, spec$step),
    spline_type = spec$type,
    spline_df = spec$df,
    spline_degree = ifelse(is.null(spec$degree), 3, spec$degree),
    spline_k = ifelse(is.null(spec$k), 10, spec$k)
  )
  
  all_results[[spec$name]] <- results
  
  # Store model fit statistics
  model_comparison <- bind_rows(
    model_comparison,
    tibble(
      Specification = spec$name,
      `R² (Conditional)` = unique(results$partial_effect$r2),
      AIC = unique(results$partial_effect$aic),
      BIC = unique(results$partial_effect$bic),
      Type = spec$type
    )
  )
}

# --- Display Model Comparison Table ---
cat("\n\n", rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON TABLE\n")
cat(rep("=", 70), "\n", sep = "")
print(model_comparison %>% arrange(desc(`R² (Conditional)`)), n = Inf)

# Prepare data for plotting - scale all curves
all_scaled_data <- map_dfr(names(all_results), function(spec_name) {
  result <- all_results[[spec_name]]
  
  result$partial_effect %>%
    mutate(
      # Calculate the sum of all predicted flows
      flow_sum = sum(pred_flow),
      # Scale prediction AND CI bounds by dividing by the sum
      pred_flow_scaled = pred_flow / flow_sum,
      ci_lower_scaled = ci_lower / flow_sum,
      ci_upper_scaled = ci_upper / flow_sum,
      specification = spec_name,
      spline_type = result$spline_type
    )
})

# Create combined plot with ALL specifications overlaid (all in navy)
p_combined <- ggplot(all_scaled_data, aes(x = distance, group = specification)) +
  geom_ribbon(
    aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
    fill = "navy",
    alpha = 0.1,
    color = NA
  ) +
  geom_line(
    aes(y = pred_flow_scaled),
    linewidth = 1,
    color = "navy",
    alpha = 0.7
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035)
  ) +
  scale_x_continuous(
    breaks = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5)
  ) +
  coord_cartesian(ylim = c(0.0015, 0.035), xlim = c(0.6, 11.4)) +
  labs(
    #title = "Robustness Check: All Spline Specifications Overlaid",
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

print(p_combined)

# Create faceted plot for detailed comparison (all in navy) - 3 columns x 2 rows
p_faceted <- ggplot(all_scaled_data, aes(x = distance)) +
  geom_ribbon(
    aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
    fill = "navy",
    alpha = 0.2
  ) +
  geom_line(
    aes(y = pred_flow_scaled),
    linewidth = 1,
    color = "navy"
  ) +
  facet_wrap(~ specification, ncol = 3, nrow = 2) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.45)
  ) +
  scale_x_continuous(
    breaks = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5)
  ) +
  coord_cartesian(ylim = c(0.0015, 0.043), xlim = c(0.6, 11.4)) +
  labs(
    #title = "Robustness Check: Individual Specifications",
    x = "Distance (km)",
    y = "Share of Flow (% of total)"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 9, face = "bold")
  )

print(p_faceted)

# Create individual plots for each specification
cat("\n\n", rep("=", 70), "\n", sep = "")
cat("GENERATING INDIVIDUAL PLOTS FOR EACH SPECIFICATION\n")
cat(rep("=", 70), "\n", sep = "")

individual_plots <- list()

for (spec_name in names(all_results)) {
  cat("\nPlotting:", spec_name, "\n")
  
  spec_data <- all_scaled_data %>%
    filter(specification == spec_name)
  
  p_individual <- ggplot(spec_data, aes(x = distance)) +
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
      breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045)
    ) +
    scale_x_continuous(
      breaks = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5)
    ) +
    coord_cartesian(ylim = c(0.0015, 0.043), xlim = c(0.6, 11.4)) +
    labs(
      title = spec_name,
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
  
  print(p_individual)
  individual_plots[[spec_name]] <- p_individual
}

# Save results
saveRDS(all_results, "~/imperial/Reach of Last Mile/dubai_robustness_all_curves.rds")
saveRDS(model_comparison, "~/imperial/Reach of Last Mile/dubai_model_comparison.rds")
saveRDS(individual_plots, "~/imperial/Reach of Last Mile/dubai_individual_plots.rds")

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("Robustness analysis complete!\n")
cat("Results saved to:\n")
cat("  - dubai_robustness_all_curves.rds\n")
cat("  - dubai_model_comparison.rds\n")
cat("  - dubai_individual_plots.rds\n")
cat(rep("=", 70), "\n", sep = "")