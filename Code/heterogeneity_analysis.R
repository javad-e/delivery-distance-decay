library(tidyverse)
library(lme4)
library(splines)
library(scales)

# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR    <- "~/imperial/Reach of Last Mile"
GROUPED_DIR <- file.path(DATA_DIR, "grouped_exports")

CONTROL_VARIABLES <- c(
  "num_vendors_dest",
  "income_cat_origin",
  "income_cat_dest",
  "vendor_entropy_origin",
  "pct_incentivized",
  "population_dest",
  "avg_delivery_fee"
)

PLOT_THEME <- theme_minimal() +
  theme(
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.major   = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor   = element_line(color = "grey90", linewidth = 0.2),
    axis.title         = element_text(size = 17),
    axis.text          = element_text(size = 16),
    legend.position    = c(1, 1),
    legend.justification = c("right", "top"),
    legend.text        = element_text(size = 15),
    legend.title       = element_text(face = "bold", size = 15)
  )

# ── Core modelling function ────────────────────────────────────────────────────

run_spline_model <- function(df,
                             control_vars = CONTROL_VARIABLES,
                             dist_min     = 0,
                             dist_max     = 18,
                             spline_step  = 2,
                             pred_points  = 200) {
  
  df <- df %>%
    filter(euclidean_km >= dist_min, euclidean_km <= dist_max) %>%
    mutate(
      origin      = factor(origin),
      destination = factor(destination),
      log_flow    = log(flow + 1)
    )
  
  available_controls <- intersect(control_vars, names(df))
  missing_controls   <- setdiff(control_vars, names(df))
  
  if (length(missing_controls) > 0)
    warning("Skipping missing controls: ", paste(missing_controls, collapse = ", "))
  
  for (var in available_controls) {
    if (is.numeric(df[[var]])) {
      df[[paste0(var, "_std")]] <- scale(df[[var]])[, 1]
    } else {
      warning(var, " is not numeric and will be skipped")
      available_controls <- setdiff(available_controls, var)
    }
  }
  
  knots        <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- ns(df$euclidean_km, knots = knots, intercept = FALSE,
                     Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  df <- bind_cols(df, as.data.frame(spline_basis))
  
  spline_terms  <- paste(colnames(spline_basis), collapse = " + ")
  control_terms <- if (length(available_controls) > 0)
    paste0(" + ", paste(paste0(available_controls, "_std"), collapse = " + "))
  else ""
  
  formula_text <- paste0("log_flow ~ 1 + ", spline_terms, control_terms,
                         " + (1 | origin) + (1 | destination)")
  cat("Formula:", formula_text, "\n\n")
  
  model <- lmer(as.formula(formula_text), data = df)
  
  distance_grid      <- tibble(distance = seq(dist_min, dist_max, length.out = pred_points))
  spline_basis_grid  <- ns(distance_grid$distance, knots = knots, intercept = FALSE,
                           Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis_grid) <- colnames(spline_basis)
  distance_grid <- bind_cols(distance_grid, as.data.frame(spline_basis_grid))
  distance_grid$origin      <- factor(levels(df$origin)[1],      levels = levels(df$origin))
  distance_grid$destination <- factor(levels(df$destination)[1], levels = levels(df$destination))
  
  for (var in available_controls)
    distance_grid[[paste0(var, "_std")]] <- 0
  
  X_grid    <- model.matrix(reformulate(c(colnames(spline_basis),
                                          paste0(available_controls, "_std"))),
                            data = distance_grid)
  beta      <- fixef(model)
  vcov_beta <- as.matrix(vcov(model))
  
  distance_grid <- distance_grid %>%
    mutate(
      pred_log_flow = as.vector(X_grid %*% beta),
      se_log_flow   = sqrt(diag(X_grid %*% vcov_beta %*% t(X_grid))),
      log_ci_lower  = pred_log_flow - 1.96 * se_log_flow,
      log_ci_upper  = pred_log_flow + 1.96 * se_log_flow,
      pred_flow     = pmax(exp(pred_log_flow) - 1, 0),
      ci_lower      = pmax(exp(log_ci_lower)  - 1, 0),
      ci_upper      = pmax(exp(log_ci_upper)  - 1, 0)
    )
  
  df$pred_full   <- predict(model, re.form = NULL)
  ss_total       <- sum((df$log_flow - mean(df$log_flow))^2)
  r2_conditional <- 1 - sum((df$log_flow - df$pred_full)^2) / ss_total
  
  list(
    partial_effect   = distance_grid %>%
      dplyr::select(distance, pred_flow, ci_lower, ci_upper) %>%
      mutate(r2 = r2_conditional),
    raw_data         = dplyr::select(df, distance = euclidean_km, flow),
    model            = model,
    control_vars_used = available_controls
  )
}

# ── Shared helpers ─────────────────────────────────────────────────────────────

scale_to_share <- function(df, group_col) {
  df %>%
    group_by(across(all_of(group_col))) %>%
    mutate(
      flow_sum          = sum(pred_flow),
      pred_flow_scaled  = pred_flow  / flow_sum,
      ci_lower_scaled   = ci_lower   / flow_sum,
      ci_upper_scaled   = ci_upper   / flow_sum
    ) %>%
    ungroup()
}

plot_share_curve <- function(df, group_col, colors, legend_title,
                             labels   = NULL,
                             ylim     = c(0.0017, 0.045),
                             xlim     = c(0.6, 12),
                             y_breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045)) {
  ggplot(df, aes(x = distance,
                 color = .data[[group_col]],
                 fill  = .data[[group_col]])) +
    geom_ribbon(aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
                alpha = 0.2, color = NA) +
    geom_line(aes(y = pred_flow_scaled), linewidth = 1.5) +
    scale_color_manual(values = colors, name = legend_title, labels = labels %||% waiver()) +
    scale_fill_manual( values = colors, name = legend_title, labels = labels %||% waiver()) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1), breaks = y_breaks) +
    scale_x_continuous(breaks = 1:12) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    labs(x = "Distance (km)", y = "Share of Flow (% of total)") +
    PLOT_THEME
}

print_r2 <- function(results_list) {
  for (nm in names(results_list)) {
    cat(nm, "— R² (conditional):",
        round(unique(results_list[[nm]]$partial_effect$r2), 4), "\n")
  }
}

# ── 1. Vendor diversity (entropy) ─────────────────────────────────────────────

cat("\n=== VENDOR DIVERSITY ===\n")

entropy_controls <- setdiff(CONTROL_VARIABLES, c("num_vendors_dest", "avg_delivery_fee",
                                                 "income_cat_dest"))

df_main   <- read_csv(file.path(DATA_DIR, "platform1_dubai_curve_estimation.csv"),
                      show_col_types = FALSE)
median_entropy <- median(df_main$vendor_entropy_destination, na.rm = TRUE)

entropy_results <- list(
  "Below median diversity" = run_spline_model(
    filter(df_main, vendor_entropy_destination <= median_entropy),
    control_vars = entropy_controls
  ),
  "Above median diversity" = run_spline_model(
    filter(df_main, vendor_entropy_destination >  median_entropy),
    control_vars = entropy_controls
  )
)

entropy_combined <- bind_rows(
  imap(entropy_results, ~ mutate(.x$partial_effect, vendor_diversity = .y))
) %>%
  mutate(vendor_diversity = factor(vendor_diversity,
                                   levels = c("Below median diversity",
                                              "Above median diversity")))

entropy_scaled <- scale_to_share(entropy_combined, "vendor_diversity")
entropy_colors <- c("Below median diversity" = "#8c510a",
                    "Above median diversity" = "#01665e")

print(plot_share_curve(entropy_scaled, "vendor_diversity", entropy_colors, "Vendor Diversity",
                       ylim = c(0.0017, 0.05),
                       y_breaks = seq(0, 0.05, 0.005)))
print_r2(entropy_results)

# ── 2. Destination income ──────────────────────────────────────────────────────

cat("\n=== DESTINATION INCOME ===\n")

income_controls <- setdiff(CONTROL_VARIABLES, "income_cat_dest")

income_results <- list(
  "Low Income"         = run_spline_model(filter(df_main, income_cat_dest == 1),
                                          control_vars = income_controls),
  "Other Income Groups" = run_spline_model(filter(df_main, income_cat_dest != 1),
                                           control_vars = income_controls)
)

income_combined <- bind_rows(
  imap(income_results, ~ mutate(.x$partial_effect, group = .y))
) %>%
  mutate(group = factor(group, levels = c("Low Income", "Other Income Groups")))

income_scaled <- scale_to_share(income_combined, "group")
income_colors <- c("Low Income" = "#2E86AB", "Other Income Groups" = "darkorchid4")

print(plot_share_curve(income_scaled, "group", income_colors, "Customer Income Group",
                       ylim = c(0.0017, 0.05),
                       y_breaks = seq(0, 0.05, 0.005)))
print_r2(income_results)
saveRDS(income_scaled, file.path(DATA_DIR, "dubai_income_curves.rds"))

# ── 3. Quarter ────────────────────────────────────────────────────────────────

cat("\n=== QUARTER ===\n")

quarter_results <- setNames(
  lapply(c("Q1", "Q2", "Q3", "Q4"), function(q) {
    cat("\n--- Processing", q, "---\n")
    df <- read_csv(file.path(GROUPED_DIR, paste0("dubai_", q, ".csv")),
                   show_col_types = FALSE)
    run_spline_model(df)
  }),
  c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec")
)

quarter_combined <- bind_rows(
  imap(quarter_results, ~ mutate(.x$partial_effect, quarter = .y))
) %>%
  mutate(quarter = factor(quarter, levels = c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec")))

quarter_scaled <- scale_to_share(quarter_combined, "quarter")
quarter_colors <- c("Jan-Mar" = "#5B7C99", "Apr-Jun" = "#7B9E7E",
                    "Jul-Sep" = "#C4976C",  "Oct-Dec" = "#8B6F7C")

print(plot_share_curve(quarter_scaled, "quarter", quarter_colors, "Time of Year",
                       ylim = c(0.0017, 0.04),
                       y_breaks = seq(0, 0.04, 0.005)))
print_r2(quarter_results)
saveRDS(quarter_scaled, file.path(DATA_DIR, "dubai_quarter_curves.rds"))

# ── 4. Package type (vertical) ────────────────────────────────────────────────

cat("\n=== PACKAGE TYPE ===\n")

vertical_results <- list(
  "Nonfood" = run_spline_model(
    read_csv(file.path(GROUPED_DIR, "dubai_workday.csv"), show_col_types = FALSE)
  ),
  "Food" = run_spline_model(
    read_csv(file.path(GROUPED_DIR, "dubai_Food.csv"),    show_col_types = FALSE)
  )
)

vertical_combined <- bind_rows(
  imap(vertical_results, ~ mutate(.x$partial_effect, day_type = .y))
) %>%
  mutate(day_type = factor(day_type, levels = c("Nonfood", "Food")))

vertical_scaled <- scale_to_share(vertical_combined, "day_type")
vertical_colors <- c("Nonfood" = "#c45002", "Food" = "#4d8d8f")

print(plot_share_curve(vertical_scaled, "day_type", vertical_colors, "Package Type",
                       ylim = c(0.0017, 0.05),
                       y_breaks = seq(0, 0.05, 0.005)))
print_r2(vertical_results)

# ── 5. Time of day ────────────────────────────────────────────────────────────

cat("\n=== TIME OF DAY ===\n")

time_periods <- c("morning", "afternoon", "evening", "night")
time_labels  <- c("Morning", "Afternoon", "Evening", "Night")

time_results <- setNames(
  lapply(time_periods, function(p) {
    cat("\n--- Processing", toupper(p), "---\n")
    df <- read_csv(file.path(GROUPED_DIR, paste0("dubai_", p, ".csv")),
                   show_col_types = FALSE)
    run_spline_model(df)
  }),
  time_labels
)

time_combined <- bind_rows(
  imap(time_results, ~ mutate(.x$partial_effect, time_period = .y))
) %>%
  mutate(time_period = factor(time_period, levels = time_labels))

time_scaled <- scale_to_share(time_combined, "time_period")
time_colors <- c("Morning" = "#8C7851", "Afternoon" = "#5B7C99",
                 "Evening" = "#8B6F7C", "Night"     = "#4A5859")

print(plot_share_curve(time_scaled, "time_period", time_colors, "Time of Day",
                       ylim = c(0.0017, 0.045),
                       y_breaks = seq(0, 0.045, 0.005)))
print_r2(time_results)
saveRDS(time_scaled, file.path(DATA_DIR, "dubai_time_curves.rds"))

# ── 6. Weekday vs weekend ─────────────────────────────────────────────────────

cat("\n=== WEEKDAY VS WEEKEND ===\n")

weekend_results <- list(
  "Weekday" = run_spline_model(
    read_csv(file.path(GROUPED_DIR, "dubai_workday.csv"), show_col_types = FALSE)
  ),
  "Weekend" = run_spline_model(
    read_csv(file.path(GROUPED_DIR, "dubai_weekend.csv"), show_col_types = FALSE)
  )
)

weekend_combined <- bind_rows(
  imap(weekend_results, ~ mutate(.x$partial_effect, day_type = .y))
) %>%
  mutate(day_type = factor(day_type, levels = c("Weekday", "Weekend")))

weekend_scaled <- scale_to_share(weekend_combined, "day_type")
weekend_colors <- c("Weekday" = "#5A8F96", "Weekend" = "#7E6AA7")

print(plot_share_curve(weekend_scaled, "day_type", weekend_colors, "Day Type",
                       ylim = c(0.0017, 0.045),
                       y_breaks = seq(0, 0.045, 0.005)))
print_r2(weekend_results)
saveRDS(weekend_scaled, file.path(DATA_DIR, "dubai_weekday_curves.rds"))

# ── 7. Incentivized vs non-incentivized ───────────────────────────────────────

cat("\n=== INCENTIVES ===\n")

incentive_controls <- setdiff(CONTROL_VARIABLES, "pct_incentivized")

incentive_results <- list(
  "noincentives" = run_spline_model(
    read_csv(file.path(GROUPED_DIR, "dubai_workday.csv"),     show_col_types = FALSE),
    control_vars = incentive_controls
  ),
  "incentivized" = run_spline_model(
    read_csv(file.path(GROUPED_DIR, "dubai_incentivized.csv"), show_col_types = FALSE),
    control_vars = incentive_controls
  )
)

incentive_combined <- bind_rows(
  imap(incentive_results, ~ mutate(.x$partial_effect, incentive_group = .y))
) %>%
  mutate(incentive_group = factor(incentive_group, levels = c("noincentives", "incentivized")))

incentive_scaled  <- scale_to_share(incentive_combined, "incentive_group")
incentive_colors  <- c("noincentives" = "#32804c", "incentivized" = "#993C91")
incentive_labels  <- c("noincentives" = "Without incentive", "incentivized" = "With incentive")

print(plot_share_curve(incentive_scaled, "incentive_group",
                       incentive_colors, "Category",
                       labels   = incentive_labels,
                       ylim     = c(0.0017, 0.04),
                       xlim     = c(0.5, 10),
                       y_breaks = seq(0, 0.04, 0.005)))
print_r2(incentive_results)