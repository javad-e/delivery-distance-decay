library(tidyverse)
library(scales)
library(patchwork)

# --- Function to fit exponential model to scaled predictions ---
fit_exponential_to_city <- function(city_data, city_name, vendor_name, platform_name) {
  # Fit exponential decay: flow = a * exp(-b * distance)
  tryCatch({
    # Try default starting values first
    model_exp <- tryCatch({
      nls(
        pred_flow_scaled ~ a * exp(-b * distance),
        data = city_data,
        start = list(a = max(city_data$pred_flow_scaled), b = 0.3),
        control = nls.control(maxiter = 1000)
      )
    }, error = function(e1) {
      # If that fails, try alternative starting values
      cat("    Retrying with alternative starting values...\n")
      tryCatch({
        nls(
          pred_flow_scaled ~ a * exp(-b * distance),
          data = city_data,
          start = list(a = max(city_data$pred_flow_scaled) * 1.1, b = 0.5),
          control = nls.control(maxiter = 2000, minFactor = 1/4096)
        )
      }, error = function(e2) {
        # Try yet another set of starting values
        cat("    Retrying with third set of starting values...\n")
        nls(
          pred_flow_scaled ~ a * exp(-b * distance),
          data = city_data,
          start = list(a = mean(city_data$pred_flow_scaled), b = 0.1),
          control = nls.control(maxiter = 3000, minFactor = 1/8192)
        )
      })
    })
    
    # Get coefficients and standard errors
    model_summary <- summary(model_exp)
    coefs <- coef(model_exp)
    ses <- model_summary$coefficients[, "Std. Error"]
    
    # Calculate 95% confidence intervals
    ci_a <- coefs["a"] + c(-1.96, 1.96) * ses["a"]
    ci_b <- coefs["b"] + c(-1.96, 1.96) * ses["b"]
    
    # Get p-values
    p_values <- model_summary$coefficients[, "Pr(>|t|)"]
    
    # Generate predictions
    distance_grid <- tibble(distance = seq(min(city_data$distance), 
                                           max(city_data$distance), 
                                           length.out = 200))
    distance_grid$pred_flow_exp <- predict(model_exp, newdata = distance_grid)
    
    # Calculate R²
    ss_total <- sum((city_data$pred_flow_scaled - mean(city_data$pred_flow_scaled))^2)
    ss_res <- sum((city_data$pred_flow_scaled - predict(model_exp))^2)
    r2_exp <- 1 - ss_res / ss_total
    
    cat("✓", city_name, "(", platform_name, "): a =", round(coefs["a"], 4), 
        "±", round(ses["a"], 4),
        ", b =", round(coefs["b"], 4), "±", round(ses["b"], 4),
        ", R² =", round(r2_exp, 4), "\n")
    
    list(
      predictions = distance_grid,
      model = model_exp,
      r2 = r2_exp,
      a = coefs["a"],
      b = coefs["b"],
      se_a = ses["a"],
      se_b = ses["b"],
      ci_a_lower = ci_a[1],
      ci_a_upper = ci_a[2],
      ci_b_lower = ci_b[1],
      ci_b_upper = ci_b[2],
      p_a = p_values["a"],
      p_b = p_values["b"],
      n_obs = nrow(city_data),
      city = city_name,
      vendor = vendor_name,
      platform = platform_name
    )
  }, error = function(e) {
    cat("✗ Error fitting", city_name, "(", platform_name, "):", e$message, "\n")
    return(NULL)
  })
}

# --- Load the saved data ---
cat("\n=== Loading Data ===\n")
platform1_data <- readRDS("~/imperial/Reach of Last Mile/platform1_curves.rds")
platform2_data <- readRDS("~/imperial/Reach of Last Mile/platform2_curves.rds")
lade_data <- readRDS("~/imperial/Reach of Last Mile/lade_curves.rds")

# --- Add vendor column to distinguish datasets ---
# platform1 is Dubai only
if (!"city" %in% names(platform1_data)) {
  platform1_data <- platform1_data %>%
    mutate(city = "Dubai")
}
platform1_data <- platform1_data %>%
  mutate(vendor = "platform1", platform = "Platform #1")

# platform2 cities
platform2_data <- platform2_data %>%
  mutate(vendor = "platform2", platform = "Platform #2")

# LaDe cities
lade_data <- lade_data %>%
  mutate(vendor = "LaDe", platform = "Platform #3")

# --- Combine the datasets ---
combined_data_all <- bind_rows(platform1_data, platform2_data, lade_data)

# Combined data without platform1 for the first plot
combined_data_no_platform1 <- bind_rows(platform2_data, lade_data)

# --- Create facet labels with platform numbers ---
combined_data_all <- combined_data_all %>%
  mutate(
    city_vendor = paste0(city, " (", platform, ")")
  )

combined_data_no_platform1 <- combined_data_no_platform1 %>%
  mutate(
    city_vendor = paste0(city, " (", platform, ")")
  )

# --- Fit exponential models to each city (ALL datasets including platform1) ---
cat("\n=== Fitting Exponential Models ===\n")
exponential_fits <- combined_data_all %>%
  group_by(city, vendor, platform, city_vendor) %>%
  group_split() %>%
  map(~ {
    city_name <- unique(.x$city)
    vendor_name <- unique(.x$vendor)
    platform_name <- unique(.x$platform)
    fit_exponential_to_city(.x, city_name, vendor_name, platform_name)
  }) %>%
  compact()  # Remove NULL entries (failed fits)

cat("\nNumber of successful fits:", length(exponential_fits), "\n")

# --- Create parameter comparison table ---
cat("\n=== Creating Parameter Table ===\n")
param_table <- map_df(exponential_fits, ~ {
  tibble(
    City = .x$city,
    Vendor = .x$vendor,
    Platform = .x$platform,
    `a (intercept)` = .x$a,
    `SE(a)` = .x$se_a,
    `CI(a) lower` = .x$ci_a_lower,
    `CI(a) upper` = .x$ci_a_upper,
    `b (decay rate)` = .x$b,
    `SE(b)` = .x$se_b,
    `CI(b) lower` = .x$ci_b_lower,
    `CI(b) upper` = .x$ci_b_upper,
    `p-value(b)` = .x$p_b,
    `R²` = .x$r2,
    `N` = .x$n_obs
  )
}) %>%
  arrange(Platform, City)

# Create display version with rounding
param_table_display <- param_table %>%
  mutate(across(c(`a (intercept)`, `CI(a) lower`, `CI(a) upper`, 
                  `b (decay rate)`, `CI(b) lower`, `CI(b) upper`, `R²`), 
                ~round(.x, 4)),
         across(c(`SE(a)`, `SE(b)`), ~round(.x, 4)),
         `p-value(b)` = formatC(`p-value(b)`, format = "e", digits = 2))

print(param_table_display)

# Export table to CSV
write_csv(param_table_display, "exponential_parameters_by_city.csv")
cat("\n✓ Parameter table saved to: exponential_parameters_by_city.csv\n")

# --- Calculate summary statistics ---
cat("\n=== Summary Statistics ===\n")
summary_stats <- param_table %>%
  summarise(
    `Mean a` = mean(`a (intercept)`),
    `SD a` = sd(`a (intercept)`),
    `Mean b` = mean(`b (decay rate)`),
    `SD b` = sd(`b (decay rate)`),
    `Mean SE(b)` = mean(`SE(b)`),
    `Mean R²` = mean(`R²`),
    `Min R²` = min(`R²`),
    `Max R²` = max(`R²`)
  )
print(summary_stats)

# Calculate coefficient of variation for decay rate
cv_b <- sd(param_table$`b (decay rate)`) / mean(param_table$`b (decay rate)`) * 100
cat("\nCoefficient of Variation for decay rate (b):", round(cv_b, 2), "%\n")

# --- Build decay_plot_data BEFORE first figure so ordering is available ---
cat("\n=== Preparing Decay Plot Data ===\n")
decay_plot_data <- param_table %>%
  mutate(
    Display_Label = paste0(City, " (", Platform, ")"),
    Display_Label = fct_reorder(Display_Label, `b (decay rate)`)
  )

# Calculate overall mean for ALL cities including Dubai
overall_mean_b <- mean(decay_plot_data$`b (decay rate)`)

# --- Prepare data for plotting (WITHOUT platform1) ---
# Get exponential fits for non-platform1 cities only
exponential_fits_no_platform1 <- exponential_fits %>%
  keep(~ .x$vendor != "platform1")

combined_data_with_exp <- combined_data_no_platform1 %>%
  left_join(
    map_df(exponential_fits_no_platform1, ~ {
      tibble(
        city = .x$city,
        vendor = .x$vendor,
        distance = .x$predictions$distance,
        pred_flow_exp = .x$predictions$pred_flow_exp
      )
    }),
    by = c("city", "vendor", "distance")
  )

# --- Align city order in first figure to match second figure (ordered by decay rate b) ---
city_order_by_b <- decay_plot_data %>%
  arrange(`b (decay rate)`) %>%              # same ordering as fct_reorder in forest plot
  filter(Platform != "Platform #1") %>%      # exclude platform1 (not shown in first figure)
  mutate(city_vendor = paste0(City, " (", Platform, ")")) %>%
  pull(city_vendor)

# Apply ordering to facet variable in first figure
combined_data_with_exp <- combined_data_with_exp %>%
  mutate(city_vendor = factor(city_vendor, levels = city_order_by_b))

# --- Create faceted plot with spline and exponential curves (NO platform1) ---
cat("\n=== Creating Comparison Plot (without platform1) ===\n")
p_comparison <- ggplot(combined_data_with_exp, aes(x = distance)) +
  # Spline confidence ribbon
  geom_ribbon(
    aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
    alpha = 0.2,
    fill = "navy"
  ) +
  # Spline line
  geom_line(
    aes(y = pred_flow_scaled, color = "Spline"),
    linewidth = 0.7,
    alpha = 0.7
  ) +
  # Exponential line
  geom_line(
    aes(y = pred_flow_exp, color = "Exponential"),
    linewidth = 0.7,
    linetype = "dashed",
    na.rm = TRUE
  ) +
  facet_wrap(~ city_vendor, ncol = 3, scales = "fixed") +
  scale_color_manual(
    name = NULL,
    values = c("Spline" = "navy", "Exponential" = "#B22222")
  ) +
  scale_y_continuous(
    breaks = c(0.01, 0.02, 0.03, 0.04),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 0.045), xlim = c(0, 7.8)) +
  labs(
    x = "Distance (km)",
    y = "Share of Flow (% of total)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 3, t = 3)),
    strip.background = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.spacing = unit(0.15, "cm"),
    legend.position = "right",
    legend.text = element_text(size = 11)
  )

print(p_comparison)

# --- Save the plot ---
ggsave("multi_city_exponential_comparison.png", 
       plot = p_comparison, 
       width = 12, 
       height = 10, 
       dpi = 300, 
       bg = "white")
cat("✓ Comparison plot saved to: multi_city_exponential_comparison.png\n")

# --- Create professional academic-style decay rate comparison plot (WITH platform1) ---
cat("\n=== Creating Professional Decay Rate Comparison Plot with CIs (including platform1) ===\n")

# Get x-axis range for proper label positioning
x_range <- range(c(decay_plot_data$`CI(b) lower`, decay_plot_data$`CI(b) upper`))

# Create the professional forest plot style with CI in front
p_decay_professional <- ggplot(decay_plot_data, 
                               aes(y = Display_Label, x = `b (decay rate)`)) +
  # Point estimates FIRST (so they appear in back)
  geom_point(
    size = 4, 
    shape = 21,
    fill = "navy",
    color = "white",
    stroke = 0.8,
    alpha = 0.8
  ) +
  # CI lines (horizontal segments) SECOND (appear in front of points)
  geom_segment(
    aes(x = `CI(b) lower`, xend = `CI(b) upper`, 
        y = Display_Label, yend = Display_Label),
    linewidth = 2, 
    color = "navy",
    alpha = 0.3
  ) +
  # Taller, narrower vertical markers at CI endpoints (left) THIRD
  geom_segment(
    aes(x = `CI(b) lower`, xend = `CI(b) lower`,
        y = as.numeric(Display_Label) - 0.25, 
        yend = as.numeric(Display_Label) + 0.25),
    linewidth = 0.5,
    color = "navy",
    alpha = 0.3
  ) +
  # Taller, narrower vertical markers at CI endpoints (right) FOURTH
  geom_segment(
    aes(x = `CI(b) upper`, xend = `CI(b) upper`,
        y = as.numeric(Display_Label) - 0.25, 
        yend = as.numeric(Display_Label) + 0.25),
    linewidth = 0.5,
    color = "navy",
    alpha = 0.3
  ) +
  # Mean reference line
  geom_vline(
    xintercept = overall_mean_b,
    linetype = "dashed", 
    color = "grey40", 
    linewidth = 0.5
  ) +
  # Annotate mean - positioned inside the plot area
  annotate("text", 
           x = overall_mean_b + diff(x_range) * 0.02,
           y = 0.5,
           label = paste0("Mean: ", sprintf("%.3f", overall_mean_b)),
           hjust = 0, 
           vjust = 0,
           size = 3.5, 
           color = "grey40", 
           fontface = "italic") +
  # Clean axis labels
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    x = expression(paste("Estimated Decay Parameter")),
    y = NULL
  ) +
  # Professional theme
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 10, color = "grey20"),
    axis.text.x = element_text(size = 10, color = "grey20"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.line.x = element_line(color = "grey40", linewidth = 0.5),
    axis.ticks.x = element_line(color = "grey40", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.caption = element_text(hjust = 0, size = 9, color = "grey40", 
                                margin = margin(t = 10)),
    plot.caption.position = "plot",
    plot.margin = margin(15, 15, 10, 15)
  )

print(p_decay_professional)

# Save professional decay rate comparison
ggsave("decay_rate_comparison_professional.png", 
       plot = p_decay_professional, 
       width = 10, 
       height = 8, 
       dpi = 400, 
       bg = "white")

# --- Check for overlapping confidence intervals ---
param_table_labeled <- param_table %>%
  mutate(label = paste0(City, " (", Platform, ")"))

comparisons <- expand_grid(
  row1 = 1:nrow(param_table_labeled),
  row2 = 1:nrow(param_table_labeled)
) %>%
  filter(row1 < row2) %>%
  mutate(
    city1 = param_table_labeled$label[row1],
    city2 = param_table_labeled$label[row2],
    b1 = param_table_labeled$`b (decay rate)`[row1],
    b2 = param_table_labeled$`b (decay rate)`[row2],
    ci1_lower = param_table_labeled$`CI(b) lower`[row1],
    ci1_upper = param_table_labeled$`CI(b) upper`[row1],
    ci2_lower = param_table_labeled$`CI(b) lower`[row2],
    ci2_upper = param_table_labeled$`CI(b) upper`[row2],
    overlap = (ci1_lower <= ci2_upper) & (ci2_lower <= ci1_upper),
    likely_different = !overlap,
    diff = abs(b1 - b2)
  ) %>%
  arrange(desc(diff))

# Show pairs that likely differ significantly
sig_different <- comparisons %>%
  filter(likely_different) %>%
  select(city1, city2, b1, b2, diff)

if (nrow(sig_different) > 0) {
  cat("\nPairs with non-overlapping 95% CIs (likely significantly different):\n")
  print(sig_different, n = Inf)
} else {
  cat("\nNo pairs found with non-overlapping confidence intervals.\n")
}

# Export comparison table
write_csv(comparisons, "decay_rate_pairwise_comparisons.csv")
cat("\n✓ Pairwise comparison table saved to: decay_rate_pairwise_comparisons.csv\n")

# --- Summary ---
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Number of cities analyzed:", nrow(param_table), "\n")
cat("Mean exponential decay rate:", round(mean(param_table$`b (decay rate)`), 4), 
    "± SE:", round(sd(param_table$`b (decay rate)`)/sqrt(nrow(param_table)), 4), "\n")
cat("Standard deviation of decay rate:", round(sd(param_table$`b (decay rate)`), 4), "\n")
cat("Mean standard error of b:", round(mean(param_table$`SE(b)`), 4), "\n")
cat("Coefficient of variation:", round(cv_b, 2), "%\n")
cat("Number of significantly different pairs:", sum(comparisons$likely_different), 
    "out of", nrow(comparisons), "total comparisons\n")
cat("\nFiles created:\n")
cat("  - exponential_parameters_by_city.csv\n")
cat("  - decay_rate_pairwise_comparisons.csv\n")
cat("  - multi_city_exponential_comparison.png (without platform1)\n")
cat("  - decay_rate_comparison_professional.png (with platform1/Dubai)\n")