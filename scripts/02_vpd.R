##########################################################
# ---- Read HH data and take averages based on day and night times ----
##########################################################

# 01_vpd_daytime.R
# This script reads half-hourly data, extracts daytime average vapour pressure deficit (VPD_F_DAY),
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
data <- read.csv("../data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/FLX_CH-Oe2_FLUXNET2015_FULLSET_HH_2004-2023_1-3.csv")
head(data)      # Preview the first few rows
colnames(data)  # Display column names

# ---- Timestamp Conversion ----
# Convert TIMESTAMP_START and TIMESTAMP_END to character format (prevents scientific notation issues)
data$TIMESTAMP_START <- as.character(data$TIMESTAMP_START)
data$TIMESTAMP_END <- as.character(data$TIMESTAMP_END)

# Extract date from TIMESTAMP_START (first 8 characters: YYYYMMDD format)
data$date <- as.Date(substr(data$TIMESTAMP_START, 1, 8), format = "%Y%m%d")

# Extract hour from TIMESTAMP_START (characters 9-10: HH format)
data$hour <- as.numeric(substr(data$TIMESTAMP_START, 9, 10))

# ---- Filter for Daytime Hours ----
# Define daytime as 6 AM to 6 PM (06:00 to 17:59)
daytime_data <- data %>%
    filter(hour >= 6 & hour < 18)

# ---- Data Aggregation ----
# Calculate daytime average VPD by date
data_daily <- daytime_data %>%
    group_by(date) %>%
    summarise(
        vpd_day = mean(VPD_F, na.rm = TRUE)  # Daytime average vapour pressure deficit
    ) %>%
    ungroup()

# Calculate nighttime average VPD (18:00 to 05:59)
nighttime_data <- data %>%
    filter(hour < 6 | hour >= 18)

nighttime_daily <- nighttime_data %>%
    group_by(date) %>%
    summarise(
        vpd_night = mean(VPD_F, na.rm = TRUE)  # Nighttime average vapour pressure deficit
    ) %>%
    ungroup()

# Calculate 24-hour average VPD (all hours)
daily_24h <- data %>%
    group_by(date) %>%
    summarise(
        vpd_24h = mean(VPD_F, na.rm = TRUE)  # 24-hour average vapour pressure deficit
    ) %>%
    ungroup()

# Merge all VPD calculations into one comprehensive dataframe
data_daily <- data_daily %>%
    left_join(nighttime_daily, by = "date") %>%
    left_join(daily_24h, by = "date")

# ---- Output ----
# Write the daily VPD data to a CSV file
write.csv(data_daily, "../data/01_data_prep/02_vpd.csv", row.names = FALSE)

# ---- Quick Checks ----
head(data_daily)    # Preview the first few rows of daily aggregated data
nrow(data_daily)    # Display total number of days in the dataset

##########################################################
# ---- Read Daily Data ----
##########################################################
# Read the daily measurement file (CSV) for comparison/validation
daily_flux_data <- read.csv("../data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/FLX_CH-Oe2_FLUXNET2015_FULLSET_DD_2004-2023_1-3.csv")

# Convert TIMESTAMP to date format for merging
daily_flux_data$TIMESTAMP <- as.character(daily_flux_data$TIMESTAMP)
daily_flux_data$date <- as.Date(daily_flux_data$TIMESTAMP, format = "%Y%m%d")

# Extract VPD_F from daily flux data for comparison with calculated values
daily_vpd_from_dd <- daily_flux_data %>%
    select(date, VPD_F) %>%
    rename(vpd_daily_flux = VPD_F)  # Rename to distinguish from calculated values

# Merge daily flux VPD with existing calculated VPD data
data_daily <- data_daily %>%
    left_join(daily_vpd_from_dd, by = "date")

# Update the CSV file with the new comparison column
write.csv(data_daily, "../data/01_data_prep/02_vpd.csv", row.names = FALSE)

# Preview the final dataset with all VPD variables
head(data_daily)




# Load ggplot2 for plotting
library(ggplot2)

# Add year column for plotting
data_daily$year <- year(data_daily$date)

# Reshape data to long format for plotting multiple variables
library(tidyr)
data_long <- data_daily %>%
    select(date, year, vpd_day, vpd_night, vpd_24h, vpd_daily_flux) %>%
    pivot_longer(cols = c(vpd_day, vpd_night, vpd_24h, vpd_daily_flux), 
                 names_to = "vpd_type", 
                 values_to = "vpd_value")

# Create a single plot with different colors for each VPD type
plot <- ggplot(data_long, aes(x = date, y = vpd_value, color = vpd_type)) +
            geom_line(alpha = 0.7, size = 0.5) +
            scale_color_manual(values = c("vpd_day" = "red", 
                                        "vpd_night" = "blue", 
                                        "vpd_24h" = "green", 
                                        "vpd_daily_flux" = "purple"),
                            labels = c("Day VPD", "Night VPD", "24h VPD", "Daily Flux VPD")) +
            scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
            labs(title = "Comparison of Different VPD Variables Over Time",
                x = "Year",
                y = "VPD (hPa)",
                color = "VPD Type") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "bottom");plot

# Create plot for Day VPD and Night VPD comparison
data_day_night <- data_daily %>%
    select(date, year, vpd_day, vpd_night) %>%
    pivot_longer(cols = c(vpd_day, vpd_night), 
                    names_to = "vpd_type", 
                    values_to = "vpd_value")

plot_day_night <- ggplot(data_day_night, aes(x = date, y = vpd_value, color = vpd_type)) +
    geom_line(alpha = 0.7, size = 0.5) +
    scale_color_manual(values = c("vpd_day" = "red", "vpd_night" = "blue"),
                        labels = c("Day VPD", "Night VPD")) +
    scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
    labs(title = "Day VPD vs Night VPD Over Time",
            x = "Year",
            y = "VPD (hPa)",
            color = "VPD Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")

# Create plot for 24h VPD and Daily Flux VPD comparison
data_24h_flux <- data_daily %>%
    select(date, year, vpd_24h, vpd_daily_flux) %>%
    pivot_longer(cols = c(vpd_24h, vpd_daily_flux), 
                    names_to = "vpd_type", 
                    values_to = "vpd_value")

plot_24h_flux <- ggplot(data_24h_flux, aes(x = date, y = vpd_value, color = vpd_type)) +
    geom_line(alpha = 0.7, size = 0.5) +
    scale_color_manual(values = c("vpd_24h" = "green", "vpd_daily_flux" = "purple"),
                        labels = c("24h VPD", "Daily Flux VPD")) +
    scale_x_date(date_labels = "%Y", date_breaks = "1 years") +
    labs(title = "24h VPD vs Daily Flux VPD Over Time",
            x = "Year",
            y = "VPD (hPa)",
            color = "VPD Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")

# Display the plots
print(plot_day_night)
print(plot_24h_flux)

# Stack all plots together for a comprehensive view
plot_combined <- gridExtra::grid.arrange(plot, plot_day_night, plot_24h_flux, ncol = 1)
plot_combined
# Create the directory for trend plots if it does not exist
if (!dir.exists("../data/01_data_prep/trend_plots")) {
    dir.create("../data/01_data_prep/trend_plots", recursive = TRUE)
}

# Save the combined plot as a PNG file
ggsave("../data/01_data_prep/trend_plots/02a_vpd_comparison_plots.png", plot = plot_combined, width = 20, height = 15, bg = "white", dpi = 300)

# stack all single plots together
library(gridExtra)
plots_stacked <- gridExtra::arrangeGrob(plot, plot_day_night, plot_24h_flux, ncol = 1)
plots_stacked
# Save the combined plot as a PNG file
ggsave("../data/01_data_prep/trend_plots/02b_vpd_combined_plots.png", plot = plots_stacked, width = 20, height = 15, bg = "white", dpi = 300)