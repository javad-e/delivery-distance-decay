# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
df_violin <- read.csv("~/imperial/Reach of Last Mile/grouped_exports/distance_violin_data.csv")
cat("Data loaded:", nrow(df_violin), "observations\n")

# Prepare labels
df_violin$incentive_label <- factor(
  df_violin$incentive,
  levels = c("no_incentive", "with_incentive"),
  labels = c("No Incentive", "Increased Distance")
)
# Shortened label - just "Quantile X"
df_violin$group_label <- paste("Quantile", df_violin$group)

# Calculate means by group and incentive
cat("Calculating means...\n")
mean_data <- df_violin %>%
  group_by(group_label, incentive_label) %>%
  summarise(
    mean_dist = mean(distance_km, na.rm = TRUE),
    se = sd(distance_km, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Calculate baseline and increase
cat("Preparing stacked bar data...\n")
stacked_data <- mean_data %>%
  group_by(group_label) %>%
  summarise(
    baseline = mean_dist[incentive_label == "No Incentive"],
    with_incentive = mean_dist[incentive_label == "Increased Distance"],
    .groups = 'drop'
  ) %>%
  mutate(
    increase = with_incentive - baseline,
    pct_change = ((with_incentive / baseline) - 1) * 100
  )

# Reshape for stacking
stacked_long <- bind_rows(
  stacked_data %>% 
    select(group_label, value = baseline) %>% 
    mutate(component = "Baseline (No Incentive)"),
  stacked_data %>% 
    select(group_label, value = increase) %>% 
    mutate(component = "Increased Distance")
) %>%
  mutate(
    component = factor(
      component,
      levels = c("Increased Distance", "Baseline (No Incentive)")
    )
  )

# Add pct_change back for labeling
stacked_long <- stacked_long %>%
  left_join(stacked_data %>% select(group_label, pct_change), by = "group_label")

print(stacked_data)

# Colors
color_baseline <- "#5DA07A"      # Green for baseline (no incentive)
color_increase <- "#BB74B4"      # Purple/pink for increase (incentive effect)

# Create stacked bar chart
cat("\nCreating stacked bar chart...\n")
p <- ggplot(stacked_long, aes(x = group_label, y = value, fill = component)) +
  
  # Stacked bars
  geom_bar(
    stat = "identity",
    width = 0.67,
    alpha = .95,
    color = "gray",
    linewidth = 0.2
  ) +
  
  # Add percentage change labels ONLY (in the middle of purple section on top)
  geom_text(
    data = stacked_data %>% filter(increase > 0),
    aes(x = group_label, y = baseline + (increase / 2),
        label = sprintf("+%.1f%%", pct_change)),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "white"
  ) +
  
  # Colors - explicit mapping
  scale_fill_manual(
    values = c(
      "Baseline (No Incentive)" = color_baseline,
      "Increased Distance" = color_increase
    ),
    name = "Category",
    breaks = c("Baseline (No Incentive)", "Increased Distance")
  ) +
  
  # Labels
  labs(
    x = "Customer Average Distance Quantile",
    y = "Distance (km)"
  ) +
  
  # Expand y-axis to fit bars properly
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  coord_cartesian(ylim = c(1.7, 6.9)) +
  
  # Theme
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 16),
    legend.position = c(1, -1),
    legend.justification = c("right", "top"),
    legend.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 15)  
  )

print(p)
