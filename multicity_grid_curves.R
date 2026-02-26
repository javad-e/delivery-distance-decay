library(tidyverse)
library(scales)

# --- Load the saved data ---
lade_data <- readRDS("~/imperial/Reach of Last Mile/platform3_curves.rds")
platform2_data <- readRDS("~/imperial/Reach of Last Mile/platform2_curves.rds")
platform1_data <- readRDS("~/imperial/Reach of Last Mile/platform1_curves.rds")

# --- Add vendor column to distinguish datasets ---
lade_data <- lade_data %>%
  mutate(vendor = "LaDe")

platform2_data <- platform2_data %>%
  mutate(vendor = "platform2")

platform1_data <- platform1_data %>%
  mutate(vendor = "platform1")

# --- Combine the datasets ---
combined_data <- bind_rows(lade_data, platform2_data, platform1_data)

# --- Create facet labels with data provider names ---
# Check if platform1 has city names
cat("platform1 cities:\n")
print(unique(platform1_data$city))

# Use factor levels to control order - Platform #1 first
combined_data <- combined_data %>%
  mutate(
    # First ensure Dubai is properly set for platform1
    city = ifelse(vendor == "platform1" & (is.na(city) | city == ""), "Dubai", city),
    city_vendor = case_when(
      vendor == "platform1" ~ paste0(city, " (Platform #1)"),
      vendor == "platform2" ~ paste0(city, " (Platform #2)"),
      vendor == "LaDe" ~ paste0(city, " (Platform #3)"),
      TRUE ~ city
    )
  ) %>%
  mutate(
    # Create factor with Platform #1 first
    city_vendor = factor(city_vendor, 
                         levels = c(
                           grep("Platform #1", unique(city_vendor), value = TRUE),
                           grep("Platform #2", unique(city_vendor), value = TRUE),
                           grep("Platform #3", unique(city_vendor), value = TRUE)
                         ))
  )

# --- Extract all unique cities from the combined dataset ---
all_cities <- combined_data %>%
  pull(city) %>%
  unique() %>%
  sort()

# Verify number of cities
cat("Number of cities found:", length(all_cities), "\n")
cat("Cities:", paste(all_cities, collapse = ", "), "\n")

# --- Create color palette for all cities ---
# Using a professional, muted color palette suitable for business/academic presentations
city_colors <- setNames(
  c("#2E4057", "#048A81", "#54457F", "#D4A574", "#8B4513",
    "#4A7C7E", "#B85042", "#6B7AA1", "#C49F47"),
  all_cities
)

# --- Calculate max predicted flow for y-axis limit (unscaled) ---
max_pred_flow <- max(combined_data$pred_flow, na.rm = TRUE)

# --- Create faceted plot in 3 rows with confidence intervals (UNSCALED) ---
p_faceted_unscaled <- ggplot(combined_data, aes(x = distance, y = pred_flow)) +
  # Add confidence ribbon (navy color)
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.2,
    fill = "navy"  # Navy color
  ) +
  # Add line (navy color)
  geom_line(
    linewidth = 0.8,
    alpha = 0.9,
    color = "navy"  # Navy color
  ) +
  # Facet by city_vendor in 3 rows with FREE y-scales
  facet_wrap(~ city_vendor, nrow = 3, scales = "free_y") +
  scale_y_continuous(
    labels = label_number(scale = 1/1000, suffix = "k", accuracy = 0.1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_x_continuous(
    breaks = c(1,2,3,4,5,6,7,8),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = c(0, 7.8)) +
  labs(
    x = "Distance (km)",
    y = "Predicted Flow (absolute count)"
  ) +
  theme_minimal(base_size = 13) +  # Smaller base size for paper
  theme(
    
    # Strip (facet label) styling - larger for city names
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 3, t = 3)),
    strip.background = element_blank(),  # Remove background color
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    
    # Grid styling
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),
    
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    
    # Minimal margins and spacing
    plot.margin = margin(3, 3, 3, 3),
    panel.spacing = unit(0.15, "cm"),
  )

# Display the faceted plot with unscaled values
print(p_faceted_unscaled)

# Save faceted plot with unscaled values
ggsave("distance_flow_by_city_unscaled.png", 
       plot = p_faceted_unscaled, 
       width = 10, 
       height = 7, 
       dpi = 300,
       bg = "white")


# --- ORIGINAL SCALED VERSION (for comparison) ---
# --- Create faceted plot in 3 rows with confidence intervals (SCALED) ---
p_faceted_scaled <- ggplot(combined_data, aes(x = distance, y = pred_flow_scaled)) +
  # Add confidence ribbon (navy color)
  geom_ribbon(
    aes(ymin = ci_lower_scaled, ymax = ci_upper_scaled),
    alpha = 0.2,
    fill = "navy"  # Navy color
  ) +
  # Add line (navy color)
  geom_line(
    linewidth = 0.8,
    alpha = 0.9,
    color = "navy"  # Navy color
  ) +
  # Facet by city_vendor in 3 rows
  facet_wrap(~ city_vendor, nrow = 3, scales = "fixed") +
  scale_y_continuous(
    breaks = c(0.01, 0.02, 0.03, 0.04),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_x_continuous(
    breaks = c(1,2,3,4,5,6,7,8),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 0.045), xlim = c(0, 7.8)) +
  labs(
    x = "Distance (km)",
    y = "Predicted Flow (% of total)"
  ) +
  theme_minimal(base_size = 13) +  # Smaller base size for paper
  theme(
    
    # Strip (facet label) styling - larger for city names
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 3, t = 3)),
    strip.background = element_blank(),  # Remove background color
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    
    # Grid styling
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.1),
    
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    
    # Minimal margins and spacing
    plot.margin = margin(3, 3, 3, 3),
    panel.spacing = unit(0.15, "cm"),
  )

# Display the scaled faceted plot
print(p_faceted_scaled)

# Save faceted plot with scaled values (original version)
ggsave("distance_flow_by_city_scaled.png", 
       plot = p_faceted_scaled, 
       width = 10, 
       height = 7, 
       dpi = 300,
       bg = "white")