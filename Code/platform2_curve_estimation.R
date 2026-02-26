library(tidyverse)
library(lme4)
library(splines)
library(scales)

# Fit spline model and return both raw data and partial effect curve
run_spline_model_return_curve <- function(df, 
                                          distance_col = "driving_km",
                                          dist_min = 0, 
                                          dist_max = 15, 
                                          spline_step = 2, 
                                          pred_points = 200) {
  df <- df %>%
    rename(distance = !!sym(distance_col)) %>%
    filter(distance >= dist_min, distance <= dist_max) %>%
    mutate(
      origin = factor(origin),
      destination = factor(destination),
      log_flow = log(flow + 1)
    )
  
  # Create spline basis
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- bs(
    df$distance,
    knots = knots,
    degree = 3,
    intercept = FALSE,
    Boundary.knots = c(dist_min, dist_max)
  )
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  
  df <- bind_cols(df, as.data.frame(spline_basis))
  
  spline_terms <- paste(paste0("spline_", seq_len(ncol(spline_basis)) - 1), 
                        collapse = " + ")
  
  formula_text <- paste0("log_flow ~ 1 + ", spline_terms,
                         " + (1 | origin) + (1 | destination)")
  
  cat("Model formula:", formula_text, "\n\n")
  
  # Fit mixed-effects model
  model <- lmer(as.formula(formula_text), data = df)
  
  # Generate predictions on grid (partial effect: fixed effects only)
  distance_grid <- tibble(distance = seq(dist_min, dist_max, length.out = pred_points))
  spline_basis_grid <- bs(
    distance_grid$distance,
    knots = knots,
    degree = 3,
    intercept = FALSE,
    Boundary.knots = c(dist_min, dist_max)
  )
  colnames(spline_basis_grid) <- paste0("spline_", seq_len(ncol(spline_basis_grid)) - 1)
  
  distance_grid <- bind_cols(distance_grid, as.data.frame(spline_basis_grid))
  distance_grid$origin <- factor(levels(df$origin)[1], levels = levels(df$origin))
  distance_grid$destination <- factor(levels(df$destination)[1], levels = levels(df$destination))
  
  # Get design matrix for fixed effects predictions
  X_grid <- model.matrix(
    reformulate(paste0("spline_", seq_len(ncol(spline_basis)) - 1)),
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
  
  list(
    partial_effect = distance_grid %>%
      dplyr::select(distance, pred_flow, ci_lower, ci_upper) %>%
      mutate(r2 = r2_conditional),
    raw_data = df %>%
      dplyr::select(distance, flow),
    model = model
  )
}

# --- CONFIGURATION: Define Your Cities Here ---
cities <- list(
  "Dubai" = "~/imperial/Reach of Last Mile/platform2_dubai_curve_estimation.csv",
  "Doha" = "~/imperial/Reach of Last Mile/platform2_doha_curve_estimation.csv",
  "Riyadh" = "~/imperial/Reach of Last Mile/platform2_riyadh_curve_estimation.csv",
  "Cairo" = "~/imperial/Reach of Last Mile/platform2_cairo_curve_estimation.csv",
  "Nairobi" = "~/imperial/Reach of Last Mile/platform2_nairobi_curve_estimation.csv",
  "Abu Dhabi" = "~/imperial/Reach of Last Mile/platform2_abudhabi_curve_estimation.csv"
)

# Choose which distance metric to use
DISTANCE_METRIC <- "h3_euclidean_km"  # Options: "driving_km", "euclidean_km", "h3_euclidean_km"

DIST_MIN <- 0
DIST_MAX <- 15  
city_colors <- c(
  "Dubai" = "#999933",
  "Doha" = "#AA4499",
  "Riyadh" = "#C39F6B",
  "Cairo" = "#8C6F72",
  "Nairobi" = "#44AA99",
  "Abu Dhabi" = "#6699CC"
)

# --- Process All Cities ---
all_results <- list()
all_partial <- list()
all_raw <- list()

for (city_name in names(cities)) {
  cat("\n=== Processing", city_name, "===\n")
  
  # Read CSV
  df <- read_csv(cities[[city_name]], show_col_types = FALSE)
  cat("Dataset dimensions:", nrow(df), "rows,", ncol(df), "columns\n")
  
  # Check if required columns exist
  required_cols <- c("origin", "destination", "flow", DISTANCE_METRIC)
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    cat("ERROR: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")
    next
  }
  
  cat("Using distance metric:", DISTANCE_METRIC, "\n")
  
  # Run spline model
  results <- run_spline_model_return_curve(
    df, 
    distance_col = DISTANCE_METRIC,
    dist_min = DIST_MIN,
    dist_max = DIST_MAX, 
    spline_step = 3.5
  )
  
  # Store results
  all_results[[city_name]] <- results
  all_partial[[city_name]] <- results$partial_effect %>% mutate(city = city_name)
  all_raw[[city_name]] <- results$raw_data %>% mutate(city = city_name)
  
  cat("Model R² (Conditional):", unique(results$partial_effect$r2), "\n")
}

# Combine all data
express_partial <- bind_rows(all_partial)
express_raw <- bind_rows(all_raw)

# --- Display R² Summary ---
cat("\n\n=== MODEL R² SUMMARY ===\n")
r2_summary <- express_partial %>%
  group_by(city) %>%
  summarise(R2 = unique(r2), .groups = "drop")
print(r2_summary)

# --- Plot 1: All Cities Combined (Original Scale) ---
p_combined <- ggplot(express_partial, aes(x = distance, color = city, fill = city)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(
    aes(y = pred_flow),
    linewidth = 1.2
  ) +
  scale_color_manual(values = city_colors) +
  scale_fill_manual(values = city_colors) +
  labs(
    title = paste("Partial Effect of Distance on Flow -", length(cities), "Cities"),
    x = "Distance (km)",
    y = "Predicted Flow",
    color = "City",
    fill = "City"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

print(p_combined)

# --- Plot 2: All Cities Combined (Scaled [0,1]) ---
express_partial_scaled_combined <- express_partial %>%
  group_by(city) %>%
  mutate(
    # Calculate the sum of all predicted flows
    flow_sum = sum(pred_flow),
    # Scale prediction AND CI bounds by dividing by the sum
    pred_flow_scaled = pred_flow / flow_sum,
    ci_lower_scaled = ci_lower / flow_sum,
    ci_upper_scaled = ci_upper / flow_sum
  ) %>%
  ungroup()

saveRDS(express_partial_scaled_combined, "~/imperial/Reach of Last Mile/platform2_curves.rds")
write.csv(express_partial_scaled_combined, "platform2_partial_scaled_combined.csv", row.names = FALSE)

# Plot (for current script)
p_combined_scaled <- ggplot(express_partial_scaled_combined, aes(x = distance, color = city, fill = city)) +
  geom_line(
    aes(y = pred_flow_scaled),
    linewidth = 0.7
  ) +
  scale_color_manual(values = city_colors) +
  scale_fill_manual(values = city_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0, max(express_partial_scaled_combined$pred_flow_scaled) * 1.1)) +
  labs(
    title = "Partial Effect of Distance on Flow - platform2 (Proportion of Total)",
    x = "Distance (km)",
    y = "Predicted Flow (proportion of total)",
    color = "City",
    fill = "City"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )
print(p_combined_scaled)

# --- Separate Plots for Each City (Original Scale) ---
for (city_name in names(cities)) {
  if (!(city_name %in% express_partial$city)) {
    cat("Skipping", city_name, "- no data available\n")
    next
  }
  
  city_data <- express_partial %>% filter(city == city_name)
  
  p_city <- ggplot(city_data, aes(x = distance)) +
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper),
      fill = city_colors[[city_name]], 
      alpha = 0.3
    ) +
    geom_line(
      aes(y = pred_flow),
      linewidth = 1.2,
      color = city_colors[[city_name]]
    ) +
    labs(
      title = paste("Partial Effect of Distance on Flow -", city_name),
      x = "Distance (km)",
      y = "Predicted Flow"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  print(p_city)
}

# --- Separate Plots for Each City (Scaled [0,1]) ---
for (city_name in names(cities)) {
  if (!(city_name %in% express_partial_scaled_combined$city)) {
    cat("Skipping", city_name, "- no data available\n")
    next
  }
  
  city_data_scaled <- express_partial_scaled_combined %>% filter(city == city_name)
  
  p_city_scaled <- ggplot(city_data_scaled, aes(x = distance)) +
    geom_ribbon(
      aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
      fill = city_colors[[city_name]], 
      alpha = 0.3
    ) +
    geom_line(
      aes(y = pred_flow_scaled),
      linewidth = 1.2,
      color = city_colors[[city_name]]
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = paste("Partial Effect of Distance on Flow -", city_name, "(Scaled [0,1])"),
      x = "Distance (km)",
      y = "Predicted Flow (scaled 0-1)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  print(p_city_scaled)
}

# --- Print detailed model summaries ---
cat("\n\n=== DETAILED MODEL SUMMARIES ===\n")
for (city_name in names(cities)) {
  if (!(city_name %in% names(all_results))) {
    next
  }
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("City:", city_name, "\n")
  cat(rep("=", 60), "\n", sep = "")
  print(summary(all_results[[city_name]]$model))
}
