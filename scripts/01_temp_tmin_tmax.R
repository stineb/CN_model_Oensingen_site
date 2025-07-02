# 01_temp_tmin_tmax.R
# This script reads half-hourly temperature data, aggregates it to daily values (mean, min, max),
# and writes the result to a CSV file.

# ---- Package Installation and Loading ----
# Ensure required packages are installed and loaded
required_packages <- c("dplyr", "lubridate")
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
    if (!pkg %in% installed_packages) {
        install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
}

# ---- Data Reading ----
# Read the half-hourly measurement file (CSV)
temp_data <- read.csv("../data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/FLX_CH-Oe2_FLUXNET2015_FULLSET_HH_2004-2023_1-3.csv")
head(temp_data)  # Preview the first few rows

# ---- Timestamp Conversion ----
# Convert TIMESTAMP_START to character (in case of scientific notation)
temp_data$TIMESTAMP_START <- as.character(temp_data$TIMESTAMP_START)
# Extract date (YYYYMMDD) and convert to Date object
temp_data$date <- as.Date(substr(temp_data$TIMESTAMP_START, 1, 8), format = "%Y%m%d")

# Extract hour from timestamp for daylight filtering
temp_data$hour <- as.numeric(substr(temp_data$TIMESTAMP_START, 9, 10))

# ---- Data Aggregation ----
# Aggregate to daily values: mean, min, and max temperature
temp_daily <- temp_data %>%
    group_by(date) %>%
    summarise(
        temp = mean(TA_F, na.rm = TRUE),  # Daily average temperature
        temp_day = mean(TA_F[hour >= 6 & hour < 18], na.rm = TRUE),  # Daily average temperature (daylight hours 6-18)
        tmin = min(TA_F, na.rm = TRUE),   # Daily minimum temperature
        tmax = max(TA_F, na.rm = TRUE)    # Daily maximum temperature
    ) %>%
    ungroup()

    # ---- Additional Data Reading for Day/Night Temperature ----
    # Read the daily measurement file for day/night temperature data
    temp_daily_file <- read.csv("../data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/FLX_CH-Oe2_FLUXNET2015_FULLSET_DD_2004-2023_1-3.csv")

    # Convert TIMESTAMP to Date format for merging
    temp_daily_file$TIMESTAMP <- as.character(temp_daily_file$TIMESTAMP)
    temp_daily_file$date <- as.Date(temp_daily_file$TIMESTAMP, format = "%Y%m%d")

    # Select relevant columns and merge with existing daily data
    temp_daily <- temp_daily %>%
        left_join(
            temp_daily_file %>% select(date, TA_F, TA_F_QC),
            by = "date"
        )


# ---- Output ----
# Write the daily aggregated data to a CSV file
write.csv(temp_daily, "../data/01_data_prep/01_temp_tmin_tmax.csv", row.names = FALSE)

# ---- Quick Checks ----
head(temp_daily)    # Preview the first few rows of daily data
nrow(temp_daily)    # Number of days in the dataset




# Load ggplot2 for plotting
library(ggplot2)
library(tidyr)

# Add year column for plotting
temp_daily$year <- year(temp_daily$date)

# Reshape data to long format for plotting daily mean, min, and max
temp_long <- temp_daily %>%
    select(date, year, temp, tmin, tmax) %>%
    pivot_longer(cols = c(temp, tmin, tmax), 
                 names_to = "temp_type", 
                 values_to = "temp_value")

# Create first plot with contrasting colors using a palette
plot1 <- ggplot(temp_long, aes(x = date, y = temp_value, color = temp_type)) +
    geom_line(alpha = 0.7, size = 0.5) +
    scale_color_manual(values = c("temp" = "#3200e6", 
                                "tmin" = "#56B4E9", 
                                "tmax" = "#D55E00"),
                      labels = c("Daily Mean", "Daily Min", "Daily Max")) +
    scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
    labs(title = "Daily Temperature: Mean, Minimum, and Maximum Over Time",
         x = "Year",
         y = "Temperature (°C)",
         color = "Temperature Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

# Create second plot for Daily Mean and Day Temperature comparison
temp_comparison <- temp_daily %>%
    select(date, temp, temp_day) %>%
    pivot_longer(cols = c(temp, temp_day), 
                 names_to = "temp_type", 
                 values_to = "temp_value")

plot2 <- ggplot(temp_comparison, aes(x = date, y = temp_value, color = temp_type)) +
    geom_line(alpha = 0.7, size = 0.5) +
    scale_color_manual(values = c("temp" = "#E69F00", 
                                "temp_day" = "#009E73"),
                      labels = c("Daily Mean (24h)", "Daily Mean (Daylight 6-18h)")) +
    scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
    labs(title = "Daily Temperature: 24h Mean vs Daylight Mean",
         x = "Year",
         y = "Temperature (°C)",
         color = "Temperature Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

# Create third plot for Daily Mean and Day Temperature comparison only
plot3 <- ggplot(temp_comparison, aes(x = date, y = temp_value, color = temp_type)) +
    geom_line(alpha = 0.7, size = 0.5) +
    scale_color_manual(values = c("temp" = "#E31A1C", 
                                "temp_day" = "#1F78B4"),
                      labels = c("Daily Mean (24h)", "Daily Mean (Daylight 6-18h)")) +
    scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
    labs(title = "Daily Temperature Comparison: 24h vs Daylight Mean",
         x = "Year",
         y = "Temperature (°C)",
         color = "Temperature Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

# Combine plots using gridExtra
library(gridExtra)
combined_plot <- grid.arrange(plot1, plot2, plot3, nrow = 3)

# Create the directory for trend plots if it does not exist
if (!dir.exists("../data/01_data_prep/trend_plots")) {
    dir.create("../data/01_data_prep/trend_plots", recursive = TRUE)
}

# Save the combined plot as a PNG file
ggsave("../data/01_data_prep/trend_plots/01a_temperature_comparison_plot.png", plot = combined_plot, width = 18, height = 10, bg = "white", dpi = 300)
