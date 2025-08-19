# Prepare drivers Data for CN Model

## ---- Package Loading ----
packages <- c("dplyr", "lubridate", "here", "tidyr", "readr", "MODISTools", "ecmwfr", "terra", "purrr")
invisible(lapply(packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
}))

## ---- Data Loading ----
### Load daily and half-hourly data once and use them for data extraction
data_dir <- here::here("data","FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3")
file_dd <- file.path(data_dir,"FLX_CH-Oe2_FLUXNET2015_FULLSET_DD_2004-2023_1-3.csv")
file_hh <- file.path(data_dir,"FLX_CH-Oe2_FLUXNET2015_FULLSET_HH_2004-2023_1-3.csv")

data_DD <- read.csv(file_dd)
data_HH <- read.csv(file_hh)

## ---- Temperature Data ----
drivers_data <- data_DD %>%
    select(TIMESTAMP, TA_F_MDS_DAY) %>%
    mutate(date = as.Date(as.character(TIMESTAMP), "%Y%m%d")) %>%
    select(date, temp = TA_F_MDS_DAY)

## ---- VPD Data (day)
### Add daily VPD from DD file for comparison
vpd_daily <- data_DD %>%
  mutate(date = as.Date(as.character(TIMESTAMP), "%Y%m%d")) %>%
  select(date, VPD_F) %>%
  mutate(vpd = VPD_F * 1.0e2) %>%
  select(date, vpd)

drivers_data <- left_join(drivers_data, vpd_daily, by = "date")


## ---- PPFD Data ----
ppfd_data <- data_DD %>%
  select(TIMESTAMP, SW_IN_F_MDS) %>%
  mutate(
    date = as.Date(as.character(TIMESTAMP), "%Y%m%d"),
    ppfd = SW_IN_F_MDS * 2.04 * 1.0e-6
  ) %>%
  select(date, ppfd)

drivers_data <- left_join(drivers_data, ppfd_data, by = "date")

## ---- Net Radiation Calculation ----

# Calculate net radiation using multiple methods and select the best available

# Check available radiation and energy balance columns
rad_cols <- grep("SW_|LW_|PPFD", names(data_DD), value = TRUE)
energy_cols <- grep("H_F|LE_F|G_F", names(data_DD), value = TRUE)

# Set constants
albedo <- 0.20
sigma <- 5.67e-8
emissivity <- 0.96

# Estimate outgoing longwave and net radiation (approximation)
data_DD <- data_DD %>%
  mutate(
    LW_OUT_est = if ("TA_F_MDS" %in% names(.)) emissivity * sigma * (TA_F_MDS + 273.15)^4 else NA_real_,
    netrad_approx = if ("TA_F_MDS" %in% names(.)) SW_IN_F * (1 - albedo) + LW_IN_F - LW_OUT_est else NA_real_,
    netrad_simple = SW_IN_F * 0.7,
    netrad_energy_balance = case_when(
      all(c("H_F_MDS", "LE_F_MDS", "G_F_MDS") %in% names(.)) ~ H_F_MDS + LE_F_MDS + G_F_MDS,
      all(c("H_F_MDS", "LE_F_MDS") %in% names(.)) ~ H_F_MDS + LE_F_MDS,
      TRUE ~ NA_real_
    )
  )

# Choose the best available method for each row
data_DD <- data_DD %>%
  mutate(
    netrad = coalesce(netrad_approx, netrad_simple, netrad_energy_balance),
    netrad_method = case_when(
      !is.na(netrad_approx) ~ "Radiation Approximation",
      !is.na(netrad_simple) ~ "Simple SW Approximation",
      !is.na(netrad_energy_balance) ~ "Energy Balance",
      TRUE ~ "No Data"
    )
  )

# Prepare net radiation dataframe for joining
netrad_df <- data_DD %>%
  mutate(date = as.Date(as.character(TIMESTAMP), "%Y%m%d")) %>%
  # select(date, netrad, netrad_method, netrad_approx, netrad_simple, netrad_energy_balance)
  select(date, netrad)

# Join net radiation to drivers_data
drivers_data <- left_join(drivers_data, netrad_df, by = "date")


# ---- Atmospheric Pressure Data ----

# Create new dataframe with date and atmospheric pressure
pa_df <- data_DD %>%
  select(TIMESTAMP, PA_F) %>%
  rename(PA_F_kpa = PA_F) %>%
  mutate(
    date = as.Date(as.character(TIMESTAMP), "%Y%m%d"),
    patm = PA_F_kpa * 1e3
  ) %>%
  # select(date, PA_F_kpa, pa)
  select(date, patm)

drivers_data <- left_join(drivers_data, pa_df, by = "date")

# ---- Precipitation Data ----
## Snow
### Precipitation as snow, in mm.
### (Is not provided explicitly in FLUXNET-type data.)
### Set snow to 0 for all dates and join to drivers_data
snow_data <- data_DD %>%
  mutate(date = as.Date(as.character(TIMESTAMP), "%Y%m%d")) %>%
  select(date) %>%
  mutate(snow = 0)

drivers_data <- left_join(drivers_data, snow_data, by = "date")

# Extract daily precipitation data and join to drivers_data
## rain

precip_data <- data_DD %>%
  select(TIMESTAMP, P_F) %>%
  mutate(
    date = as.Date(as.character(TIMESTAMP), "%Y%m%d"),
    rain = P_F / (60 * 60 * 24)  # precipitation rate (mm/s) to mm/day
  ) %>%
  select(date, rain)

drivers_data <- left_join(drivers_data, precip_data, by = "date")

# ---- Temperature Min/Max Data ----
temp_minmax <- data_HH %>%
    mutate(date = as.Date(substr(as.character(TIMESTAMP_START), 1, 8), "%Y%m%d")) %>%
    group_by(date) %>%
    summarise(
        tmin = min(TA_F, na.rm = TRUE),
        tmax = max(TA_F, na.rm = TRUE)
    ) %>%
    ungroup()

drivers_data <- left_join(drivers_data, temp_minmax, by = "date")

## ---- Wind speed Data ----
wind_data <- data_DD %>%
    select(TIMESTAMP, WS_F) %>%
    mutate(date = as.Date(as.character(TIMESTAMP), "%Y%m%d")) %>%
    select(date, vwind = WS_F)
drivers_data <- left_join(drivers_data, wind_data, by = "date")

# ----------------------------# ----------------------------# ----------------------------
#----- fapar and lai data using MODISTools------------Comment out lines below (it will take 2-3 hours to download)
# ----------------------------# ----------------------------# ----------------------------
# # ----------------------------
# # Parameters
# # ----------------------------
# lat  <- 47.28631
# lon  <- 7.734333
# start_date <- "2004-01-01"
# end_date   <- "2023-12-31"
# 
# out_dir <- "data/01_data_prep_modistools"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ----------------------------
# # Helpers
# # ----------------------------
# 
# download_modis <- function(product, bands, site_name, km_lr = 1, km_ab = 1) {
#   message(sprintf("== Downloading %s for %s ==", product, site_name))
#   dat <- MODISTools::mt_subset(
#     product  = product,
#     lat      = lat,
#     lon      = lon,
#     start    = start_date,
#     end      = end_date,
#     band     = bands,
#     km_lr    = km_lr,
#     km_ab    = km_ab,
#     internal = TRUE,
#     progress = FALSE
#   )
#   if (nrow(dat) == 0) stop(sprintf("No data returned for %s", product))
#   dat <- tibble::as_tibble(dat)
#   dat
# }
# 
# to_daily <- function(df, date_col = "calendar_date", value_col = "value") {
#   # Keep one value per date (if multiple bands, filter beforehand)
#   df <- df %>%
#     mutate(date = as.Date(.data[[date_col]])) %>%
#     arrange(date) %>%
#     distinct(date, .keep_all = TRUE)
#   
#   full_seq <- tibble(date = seq(as.Date(start_date), as.Date(end_date), by = "day"))
#   
#   full_seq %>%
#     left_join(df %>% select(date, !!value_col), by = "date") %>%
#     fill(all_of(value_col), .direction = "down") %>%
#     fill(all_of(value_col), .direction = "up") %>%
#     filter(!is.na(.data[[value_col]]))
# }
# 
# # ----------------------------
# # FAPAR: MOD15A2H
# # Bands: Fpar_500m, FparStdDev_500m, FparLai_QC
# # Fpar_500m is 0–100; convert to fraction 0–1 by dividing by 100
# # ----------------------------
# 
# fapar_raw <- download_modis(
#   product = "MOD15A2H",
#   bands   = c("Fpar_500m", "FparStdDev_500m", "FparLai_QC"),
#   site_name = "site_fapar",
#   km_lr = 1, km_ab = 1
# )
# 
# # Save the raw subset
# fapar_raw_file <- file.path(out_dir, "06a_fapar_2004-2023_raw.csv")
# write_csv(fapar_raw, fapar_raw_file)
# message(sprintf("Saved raw FAPAR subset: %s", fapar_raw_file))
# 
# # Extract main band and convert to fraction
# fapar_main <- fapar_raw %>%
#   filter(band == "Fpar_500m") %>%
#   mutate(
#     date  = as.Date(calendar_date),
#     fapar = value / 100
#   ) %>%
#   select(date, fapar) %>%
#   arrange(date)
# 
# fapar_main_file <- file.path(out_dir, "06a1_fapar_2004-2023.csv")
# write_csv(fapar_main, fapar_main_file)
# message(sprintf("Saved FAPAR (8-day) in fraction units: %s", fapar_main_file))
# 
# # Convert 8-day FAPAR to daily by LOCF + BOCF
# fapar_daily <- to_daily(df = fapar_main, date_col = "date", value_col = "fapar")
# fapar_daily_file <- file.path(out_dir, "06a2_daily_fapar_2004-2023.csv")
# write_csv(fapar_daily, fapar_daily_file)
# message(sprintf("Saved daily FAPAR: %s", fapar_daily_file))
# message(sprintf("Daily FAPAR has %d rows and %d columns", nrow(fapar_daily), ncol(fapar_daily)))
# 
# # ----------------------------
# # LAI: MOD15A2H
# # Band: Lai_500m and FparLai_QC retained for possible QA
# # Lai_500m uses scale factor 0.1 (m2/m2)
# # ----------------------------
# 
# lai_raw <- download_modis(
#   product = "MOD15A2H",
#   bands   = c("Lai_500m", "FparLai_QC"),
#   site_name = "site_lai",
#   km_lr = 1, km_ab = 1
# )
# 
# lai_raw_file <- file.path(out_dir, "06b_lai_2004-2023_raw.csv")
# write_csv(lai_raw, lai_raw_file)
# message(sprintf("Saved raw LAI subset: %s", lai_raw_file))
# 
# # Main LAI band to physical units
# lai_main <- lai_raw %>%
#   filter(band == "Lai_500m") %>%
#   mutate(
#     date = as.Date(calendar_date),
#     lai  = value * 0.1
#   ) %>%
#   select(date, lai) %>%
#   arrange(date)
# 
# lai_main_file <- file.path(out_dir, "06b1_lai_2004-2023.csv")
# write_csv(lai_main, lai_main_file)
# message(sprintf("Saved LAI (8-day) in m^2/m^2: %s", lai_main_file))
# 
# # Daily LAI with the same expansion
# lai_daily <- to_daily(df = lai_main, date_col = "date", value_col = "lai")
# lai_daily_file <- file.path(out_dir, "06b2_daily_lai_2004-2023.csv")
# write_csv(lai_daily, lai_daily_file)
# message(sprintf("Saved daily LAI: %s", lai_daily_file))
# message(sprintf("Daily LAI has %d rows and %d columns", nrow(lai_daily), ncol(lai_daily)))
# ---------------------------- read the fapar and lai data and join to drivers_data
# Read the FAPAR and LAI daily data
fapar_daily <- read_csv("data/01_data_prep_modistools/06a2_daily_fapar_2004-2023.csv")
lai_daily <- read_csv("data/01_data_prep_modistools/06b2_daily_lai_2004-2023.csv")
# left join to the drivers data
drivers_data <- left_join(drivers_data, fapar_daily, by = "date")
drivers_data <- left_join(drivers_data, lai_daily, by = "date")



# ---- CO2 Data (Daily) ----

# Extract daily CO2 data from DD file and join to drivers_data
co2_data <- data_DD %>%
    select(TIMESTAMP, CO2_F_MDS) %>%
    mutate(date = as.Date(as.character(TIMESTAMP), "%Y%m%d")) %>%
    select(date, co2 = CO2_F_MDS)

drivers_data <- left_join(drivers_data, co2_data, by = "date")


# ----- cloud coverage (ccov) data ---- Comment out lines below (it will take 12-14 hours to download)
# library(ecmwfr)
# 
# cds_user <- "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6"  # your CDS login (email or UID)
# 
# wf_set_key(user = "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6",
#            key  = "01783a93-8c21-4104-a0c5-180da2701296")
# # ============================================
# # ERA5 Cloud Cover (total_cloud_cover) Downloader & Processor
# # Requires: ecmwfr, terra, dplyr, tidyr, readr, lubridate, purrr
# # ============================================
# 
# # install.packages(c("ecmwfr","terra","dplyr","tidyr","readr","lubridate","purrr"))  # once
# library(ecmwfr)
# library(terra)
# library(dplyr)
# library(tidyr)
# library(readr)
# library(lubridate)
# library(purrr)
# 
# # --- Network hardening ---
# options(timeout = 24 * 3600)  # allow long operations overall
# # Raise libcurl low-speed window from 600s to 3600s (1 hour)
# httr::set_config(httr::config(low_speed_time = 3600, low_speed_limit = 1))
# 
# 
# # ------------------------------------------------------------
# # 0) SETTINGS
# # ------------------------------------------------------------
# # Store your CDS key once per machine:
# wf_set_key(user = "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6", key = "01783a93-8c21-4104-a0c5-180da2701296")
# cds_user    <- "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6"  # your CDS login (email or UID)
# 
# out_nc_dir  <- "data/02_data_ccov/data-raw/cloud_cover"  # NetCDF output
# out_csv_dir <- "data/02_data_ccov/cloud_cover"           # CSV output
# # dir.create(out_nc_dir,  showWarnings = FALSE, recursive = TRUE)
# # dir.create(out_csv_dir, showWarnings = FALSE, recursive = TRUE)
# 
# sites <- tibble::tibble(
#   sitename   = "Oensingen",
#   lat        = 47.3509,
#   lon        = 7.734333,
#   year_start = 2004
# )
# 
# end_year <- 2023   # last full year you want
# 
# # ----------------------------
# # Helpers
# # ----------------------------
# pad2 <- function(x) sprintf("%02d", x)
# 
# # Exact day list for a given year-month (avoids invalid days like 30/31 in short months)
# days_of_month <- function(year, month) {
#   d1 <- as.Date(sprintf("%04d-%02d-01", year, month))
#   d2 <- ceiling_date(d1, "month") - days(1)
#   pad2(seq_len(day(d2)))
# }
# 
# # One MONTH request payload with valid day list
# month_request <- function(lat, lon, target, year, month) {
#   list(
#     product_type       = "reanalysis",
#     format             = "netcdf",
#     variable           = "total_cloud_cover",
#     year               = as.character(year),
#     month              = pad2(month),
#     day                = days_of_month(year, month),
#     time               = sprintf("%02d:00", 0:23),
#     # area order: North, West, South, East (tight box around point)
#     area               = c(lat + 0.05, lon - 0.05, lat - 0.05, lon + 0.05),
#     dataset_short_name = "reanalysis-era5-single-levels",
#     target             = target
#   )
# }
# 
# # Robust single-month submit with retries and backoff
# download_one <- function(user, request, path, tries = 6, base_sleep = 15) {
#   for (i in seq_len(tries)) {
#     err <- tryCatch({
#       ecmwfr::wf_request(user = user, request = request, path = path, time_out = 6 * 3600)
#       NULL
#     }, error = function(e) conditionMessage(e))
#     
#     if (is.null(err)) return(invisible(TRUE))
#     
#     if (grepl("Timeout was reached", err, fixed = TRUE)) {
#       # back off and retry on low-speed/queue timeouts
#       Sys.sleep(base_sleep * (2^(i - 1)))
#     } else {
#       stop("CDS error for ", request$target, ": ", err)
#     }
#     
#     if (i == tries) stop("CDS error for ", request$target, ": ", err)
#   }
# }
# 
# verify_nc <- function(filepath) {
#   ok <- try({
#     r <- terra::rast(filepath)
#     nlyr(r) > 0
#   }, silent = TRUE)
#   isTRUE(ok)
# }
# 
# # Build monthly job list (skips already present files)
# build_jobs <- function(sites_tbl, end_year, out_nc_dir) {
#   jobs <- list()
#   for (i in seq_len(nrow(sites_tbl))) {
#     site <- sites_tbl$sitename[i]
#     lat  <- sites_tbl$lat[i]
#     lon  <- sites_tbl$lon[i]
#     y0   <- sites_tbl$year_start[i]
#     for (y in y0:end_year) {
#       for (m in 1:12) {
#         target <- sprintf("%s_%d_%02d.nc", site, y, m)
#         jobs[[length(jobs) + 1]] <- list(
#           site = site,
#           target = target,
#           req = month_request(lat, lon, target, y, m)
#         )
#       }
#     }
#   }
#   present <- list.files(out_nc_dir, pattern = "\\.nc$", full.names = FALSE)
#   Filter(function(j) !(j$target %in% present), jobs)
# }
# 
# # ----------------------------
# # Submit monthly jobs sequentially
# # ----------------------------
# if (is.null(ecmwfr::wf_get_key(user = cds_user))) {
#   stop("No CDS key found for '", cds_user, "'. Call wf_set_key(user = ..., key = ...) once, then retry.")
# }
# 
# jobs <- build_jobs(sites, end_year, out_nc_dir)
# 
# if (length(jobs) == 0) {
#   message("All monthly NetCDF files already present.")
# } else {
#   message("Submitting ", length(jobs), " monthly jobs (sequential).")
#   for (k in seq_along(jobs)) {
#     j <- jobs[[k]]
#     message(sprintf("[%d/%d] %s", k, length(jobs), j$target))
#     download_one(user = cds_user, request = j$req, path = out_nc_dir, tries = 6, base_sleep = 20)
#     
#     fp <- file.path(out_nc_dir, j$target)
#     if (!file.exists(fp) || !verify_nc(fp)) {
#       stop("Downloaded file failed verification: ", fp)
#     }
#     
#     # polite pause to reduce queue pressure
#     Sys.sleep(8)
#   }
# }
# 
# # ----------------------------
# # Process to daily means (per site) — robust version
# # ----------------------------
# process_site <- function(site,
#                          drop_leap = FALSE,
#                          allowed_ext = c("nc", "grib", "grb")) {
#   
#   # Match monthly files like: Oensingen_2014_01.nc / .grib / .grb
#   ext_re <- paste(allowed_ext, collapse = "|")
#   pat <- paste0("^", site, "_\\d{4}_\\d{2}\\.(", ext_re, ")$")
#   files <- list.files(out_nc_dir, pattern = pat, full.names = TRUE)
#   if (length(files) == 0) return(NULL)
#   
#   # Sort by YYYYMM from filename
#   key <- as.integer(gsub("^.*_(\\d{4})_(\\d{2})\\..*$", "\\1\\2", basename(files)))
#   files <- files[order(key)]
#   
#   out <- vector("list", length(files))
#   
#   for (i in seq_along(files)) {
#     f  <- files[i]
#     bn <- basename(f)
#     yy <- as.integer(sub("^.*_(\\d{4})_(\\d{2})\\..*$", "\\1", bn))
#     mm <- as.integer(sub("^.*_(\\d{4})_(\\d{2})\\..*$", "\\2", bn))
#     
#     # Open monthly stack of hourly total_cloud_cover
#     r  <- suppressWarnings(terra::rast(f))
#     nl <- terra::nlyr(r)
#     if (nl == 0) stop("Empty raster: ", f)
#     
#     # Expected hours in month
#     d1 <- as.Date(sprintf("%04d-%02d-01", yy, mm))
#     d2 <- lubridate::ceiling_date(d1, "month") - lubridate::days(1)
#     nd <- as.integer(lubridate::day(d2))
#     expected <- nd * 24L
#     
#     # If CDS delivered fewer/more layers, trim to full 24-hour days
#     if (nl < 24L) stop("Too few layers (<24) in ", f, " (have ", nl, ").")
#     if (nl %% 24L != 0L) {
#       nd_trim <- floor(nl / 24L)
#       warning(sprintf("Layer count mismatch for %s: have %d, expected %d. Using first %d days.",
#                       bn, nl, expected, nd_trim))
#       r  <- r[[1:(nd_trim * 24L)]]
#       nd <- nd_trim
#       d2 <- d1 + (nd - 1L)
#     }
#     
#     # Group 24-hour chunks per day and average
#     idx   <- rep(seq_len(nd), each = 24L)
#     r_day <- terra::tapp(r, index = idx, fun = "mean")
#     
#     vals  <- terra::values(r_day, dataframe = TRUE)
#     dates <- seq(d1, d2, by = "1 day")
#     ccov  <- as.numeric(vals[[1]])  # single variable
#     
#     out[[i]] <- tibble::tibble(
#       sitename = site,
#       date     = dates,
#       ccov     = ccov
#     )
#     
#     # Free memory early
#     rm(r, r_day, vals); gc()
#   }
#   
#   df <- dplyr::bind_rows(out) |> dplyr::arrange(date)
#   if (drop_leap) {
#     df <- df |>
#       dplyr::filter(!(lubridate::month(date) == 2 & lubridate::mday(date) == 29))
#   }
#   df
# }
# 
# # Run processing
# dl <- lapply(sites$sitename, process_site)
# dl <- dl[!vapply(dl, is.null, logical(1))]
# if (length(dl) == 0) stop("No processed cloud cover data produced.")
# 
# # Save per-site and combined CSVs
# invisible(lapply(dl, function(df_site) {
#   out_site <- file.path(out_csv_dir, paste0(unique(df_site$sitename), "_ccov_daily.csv"))
#   readr::write_csv(df_site, out_site)
#   message("Saved: ", out_site, " (", nrow(df_site), " rows)")
# }))
# 
# ccov_all <- dplyr::bind_rows(dl) |> dplyr::arrange(sitename, date)
# out_all  <- file.path(out_csv_dir, "ccov_daily_all_sites.csv")
# readr::write_csv(ccov_all, out_all)
# message("Saved combined: ", out_all, " (", nrow(ccov_all), " rows)")




# read the data csv
data_ccov <- read_csv("data/02_data_ccov/cloud_cover/Oensingen_ccov_daily.csv")
# select only date and ccocv columns
data_ccov <- data_ccov %>%
    select(date, ccov) 
# join the cloud coverage data to drivers_data
drivers_data <- left_join(drivers_data, data_ccov, by = "date") 


# ---- GPP, LE, and NEE Data ----

gpp_data <- data_DD %>%
    select(TIMESTAMP, GPP_NT_VUT_REF, 
           # GPP_NT_VUT_REF_QC, 
           LE_F_MDS, 
        #    LE_F_MDS_QC, 
           NEE_VUT_REF, 
           NEE_VUT_REF_QC) %>%
    mutate(date = as.Date(as.character(TIMESTAMP), format = "%Y%m%d")) %>%
    select(date, gpp = GPP_NT_VUT_REF, 
           # gpp_qc = GPP_NT_VUT_REF_QC, 
           le = LE_F_MDS, 
        #    le_qc = LE_F_MDS_QC, 
           nee = NEE_VUT_REF, 
           nee_qc = NEE_VUT_REF_QC)

drivers_data <- left_join(drivers_data, gpp_data, by = "date")
# ---- Save the drivers data ----
# Save the prepared drivers data to a CSV file
drivers_data_file <- file.path(data_dir, "drivers_data_CH-Oe2_19082025.csv")
write_csv(drivers_data, drivers_data_file)

# read rds forcing data file previously uploaded for cn model as drivers data
drivers <- readRDS("data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/rsofun_driver_data_v3.4.2.rds") %>% 
  dplyr::filter(sitename == "CH-Oe2") 

# replace forcing data forcing_df <- drivers$forcing[[1]] with a new csv file

# read csv file with forcing data
forcing_df <- read_csv("data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/drivers_data_CH-Oe2_19082025.csv") 

# remove drivers$forcing[[1]] from drivers_data
drivers <- drivers %>% 
  dplyr::select(-forcing)

# add forcing_df to drivers_data as forcing
drivers <- drivers %>% 
  dplyr::mutate(forcing = list(forcing_df))


# Add management data 
# add management data into drivers$forcing
# Management data
#N input
n_input <- read.csv("data/mangement_data/04_23choe2_ninput.csv")

#fharv 
fharv <- read.csv("data/mangement_data/04_23ch0e2_fharv.csv")
# convert to fraction from percentage to be consistent with the rest of the data
fharv$fharv <- fharv$fharv / 100

#seed 
seed <- read.csv("data/mangement_data/04_23ch0e2_seed.csv")

#yieldcn  
yieldcn <- read.csv("data/mangement_data/04_23ch0e2_yieldcn.csv")

## N deposition------------
# Reactive N input needs to be in gN per day
# added to forcing time series: specify quantity of N added on which day

# Example N input given a constant rate each day
n_input_add <- function(data, n_input) {
  # Ensure date columns are in Date format
  n_input <- n_input %>% mutate(date = as.Date(date))
  
  # Remove duplicates from n_input based on the date column
  n_input <- n_input %>% distinct(date, .keep_all = TRUE)
  
  data <- data %>%
    mutate(forcing = purrr::map(forcing, ~ {
      # Merge forcing with n_input on date
      merged <- left_join(., n_input, by = "date")
      
      # Fill NA values in dno3 and dnh4 with the default dnh4 value from the n_input DataFrame
      merged <- merged %>%
        mutate(
          dno3 = ifelse(is.na(dno3), 0.002003263, dno3),
          dnh4 = ifelse(is.na(dnh4), 0.002017981, dnh4)
        )
      
      # Ensure the merged DataFrame does not have more rows than the original forcing DataFrame
      merged <- merged %>% slice(1:nrow(.))
      
      return(merged)
    }))
  
  return(data)
}
# Add N input to forcing data this will join the drivers data with the N input data
drivers <- n_input_add(drivers, n_input)


## Harvesting-----------
# The fraction of biomass harvested per day needs to be specified in the forcing time series
# cseed and nseed new seeds added after harvesting
# Example driver update assumes harvesting is 0 and new seeds planted after harvest

fharv_seed_add <- function(data, fharv, seed) {
  # Ensure date columns are in Date format
  fharv <- fharv %>% mutate(date = as.Date(date)) %>% distinct(date, .keep_all = TRUE)
  seed <- seed %>% mutate(date = as.Date(date)) %>% distinct(date, .keep_all = TRUE)
  
  data <- data %>%
    mutate(forcing = purrr::map(forcing, ~ {
      # Merge forcing with fharv and seed on date
      merged <- left_join(., fharv, by = "date") %>%
        left_join(seed, by = "date") %>%
        mutate(
          fharv = coalesce(fharv, 0),    # Fill NA values in fharv
          cseed = coalesce(cseed, 0),    # Fill NA values in cseed
          nseed = coalesce(nseed, 0)      # Fill NA values in nseed
        ) %>%
        slice(1:nrow(.))  # Ensure merged DataFrame does not have more rows than original
      
      return(merged)
    }))
  
  return(data)
}
# Add harvesting and seed data to forcing data 
drivers <- fharv_seed_add(drivers, fharv, seed)

## Simulation parameters------------
# The spinup of cn_model must be long enough to equilibrate fluxes

# Function to modify specific columns in each dataframe
modify_params <- function(df_list, spinupyears_val, recycle_val) {
  # Map over each dataframe in the list
  df_list <- map(df_list, ~ {
    # Modify specific columns
    mutate(.x,
           spinupyears = spinupyears_val,
           recycle = recycle_val)
  })
  
  return(df_list)
}

# FLUXNET data for site CH-Oe2 (Oensingen, Switzerland)
drivers$params_siml <- modify_params(drivers$params_siml, 10, 2)
drivers$params_siml[[1]]$c_only <- TRUE

nrow(drivers$forcing[[1]])

# remove rows from drivers$forcing[[1]], when it is leap year such as 2004-02-29
drivers$forcing[[1]] <- drivers$forcing[[1]][!grepl("-02-29", drivers$forcing[[1]]$date), ]
# show number of rows in drivers$forcing[[1]]
nrow(drivers$forcing[[1]])

# save the updated drivers data
saveRDS(drivers, "data/CH-Oe2_2004-2023_final_ready_for_CNmodel_run_19082025.rds")




### Comparison of Drivers Data with _fdk data
# Load packages without warnings and messages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(visdat)
  library(here)
  library(lubridate)
  library(readr)
  library(naniar)
  library(purrr)
})

# load drivers data
drivers <- readRDS(here::here("data/CH-Oe2_2004-2023_final_ready_for_CNmodel_run_19082025.rds"))
colnames(drivers$forcing[[1]])
str(drivers$forcing[[1]])

# check drivers consistency with FluxDataKit
drivers_fdk <- read_rds("data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/rsofun_driver_data_v3.4.2.rds")
drivers_fdk <- drivers_fdk |> filter(sitename == "CH-Oe2")
colnames(drivers_fdk$forcing[[1]])
str(drivers_fdk$forcing[[1]])

# Create directory if it doesn't exist
if (!dir.exists("comparison_plots")) dir.create("comparison_plots")

# Variables to plot
variables <- c('temp', 'vpd', 'ppfd', 'patm', 'snow', 'rain', 
         'tmin', 'tmax', 'vwind', 'fapar', 'co2', 'ccov', 'nee')

# Join datasets
tmp <- drivers_fdk$forcing[[1]] |> 
  right_join(
  drivers$forcing[[1]],
  by = join_by(date),
  suffix = c("_fdk", "_muh")
  )

# Function to create and save scatter plot
create_comparison_plot <- function(var_name, data, index) {
  # Format index for filename (ensure 2 digits)
  idx <- sprintf("%02d", index)
  
  # Create scatter plot
  p1 <- data |> 
  ggplot(aes_string(paste0(var_name, "_fdk"), paste0(var_name, "_muh"))) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = paste("Comparison of", var_name),
     x = paste(var_name, "FDK"),
     y = paste(var_name, "MUH"))
  
  # Create time series plot
  p2 <- data |> 
  ggplot() +
  geom_line(aes(date, !!sym(paste0(var_name, "_muh"))), color = "red") +
  geom_line(aes(date, !!sym(paste0(var_name, "_fdk")))) +
  labs(title = paste("Time series of", var_name),
     x = "Date", 
     y = var_name)
  
  # Combine plots
  combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  
  # Save plot
  ggsave(paste0("comparison_plots/", idx, "_new_", var_name, ".png"), 
     combined_plot, 
     width = 10, height = 5)
  
  return(combined_plot)
}

# Create plots for all variables
plots <- map2(variables, 1:length(variables), 
       ~create_comparison_plot(.x, tmp, .y))

# Print message
cat("All comparison plots saved in 'comparison_plots' directory\n")
