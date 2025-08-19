library(ecmwfr)

cds_user <- "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6"  # your CDS login (email or UID)

wf_set_key(user = "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6",
           key  = "01783a93-8c21-4104-a0c5-180da2701296")
# ============================================
# ERA5 Cloud Cover (total_cloud_cover) Downloader & Processor
# Requires: ecmwfr, terra, dplyr, tidyr, readr, lubridate, purrr
# ============================================

# install.packages(c("ecmwfr","terra","dplyr","tidyr","readr","lubridate","purrr"))  # once
library(ecmwfr)
library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(purrr)

# --- Network hardening ---
options(timeout = 24 * 3600)  # allow long operations overall
# Raise libcurl low-speed window from 600s to 3600s (1 hour)
httr::set_config(httr::config(low_speed_time = 3600, low_speed_limit = 1))


# ------------------------------------------------------------
# 0) SETTINGS
# ------------------------------------------------------------
# Store your CDS key once per machine:
wf_set_key(user = "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6", key = "01783a93-8c21-4104-a0c5-180da2701296")
cds_user    <- "b9d56fdb-0653-46c4-a8fa-a0ce9188c2b6"  # your CDS login (email or UID)

out_nc_dir  <- "data/02_data_ccov/data-raw/cloud_cover"  # NetCDF output
out_csv_dir <- "data/02_data_ccov/cloud_cover"           # CSV output
# dir.create(out_nc_dir,  showWarnings = FALSE, recursive = TRUE)
# dir.create(out_csv_dir, showWarnings = FALSE, recursive = TRUE)

sites <- tibble::tibble(
  sitename   = "Oensingen",
  lat        = 47.3509,
  lon        = 7.734333,
  year_start = 2004
)

end_year <- 2023   # last full year you want

# ----------------------------
# Helpers
# ----------------------------
pad2 <- function(x) sprintf("%02d", x)

# Exact day list for a given year-month (avoids invalid days like 30/31 in short months)
days_of_month <- function(year, month) {
  d1 <- as.Date(sprintf("%04d-%02d-01", year, month))
  d2 <- ceiling_date(d1, "month") - days(1)
  pad2(seq_len(day(d2)))
}

# One MONTH request payload with valid day list
month_request <- function(lat, lon, target, year, month) {
  list(
    product_type       = "reanalysis",
    format             = "netcdf",
    variable           = "total_cloud_cover",
    year               = as.character(year),
    month              = pad2(month),
    day                = days_of_month(year, month),
    time               = sprintf("%02d:00", 0:23),
    # area order: North, West, South, East (tight box around point)
    area               = c(lat + 0.05, lon - 0.05, lat - 0.05, lon + 0.05),
    dataset_short_name = "reanalysis-era5-single-levels",
    target             = target
  )
}

# Robust single-month submit with retries and backoff
download_one <- function(user, request, path, tries = 6, base_sleep = 15) {
  for (i in seq_len(tries)) {
    err <- tryCatch({
      ecmwfr::wf_request(user = user, request = request, path = path, time_out = 6 * 3600)
      NULL
    }, error = function(e) conditionMessage(e))
    
    if (is.null(err)) return(invisible(TRUE))
    
    if (grepl("Timeout was reached", err, fixed = TRUE)) {
      # back off and retry on low-speed/queue timeouts
      Sys.sleep(base_sleep * (2^(i - 1)))
    } else {
      stop("CDS error for ", request$target, ": ", err)
    }
    
    if (i == tries) stop("CDS error for ", request$target, ": ", err)
  }
}

verify_nc <- function(filepath) {
  ok <- try({
    r <- terra::rast(filepath)
    nlyr(r) > 0
  }, silent = TRUE)
  isTRUE(ok)
}

# Build monthly job list (skips already present files)
build_jobs <- function(sites_tbl, end_year, out_nc_dir) {
  jobs <- list()
  for (i in seq_len(nrow(sites_tbl))) {
    site <- sites_tbl$sitename[i]
    lat  <- sites_tbl$lat[i]
    lon  <- sites_tbl$lon[i]
    y0   <- sites_tbl$year_start[i]
    for (y in y0:end_year) {
      for (m in 1:12) {
        target <- sprintf("%s_%d_%02d.nc", site, y, m)
        jobs[[length(jobs) + 1]] <- list(
          site = site,
          target = target,
          req = month_request(lat, lon, target, y, m)
        )
      }
    }
  }
  present <- list.files(out_nc_dir, pattern = "\\.nc$", full.names = FALSE)
  Filter(function(j) !(j$target %in% present), jobs)
}

# ----------------------------
# Submit monthly jobs sequentially
# ----------------------------
if (is.null(ecmwfr::wf_get_key(user = cds_user))) {
  stop("No CDS key found for '", cds_user, "'. Call wf_set_key(user = ..., key = ...) once, then retry.")
}

jobs <- build_jobs(sites, end_year, out_nc_dir)

if (length(jobs) == 0) {
  message("All monthly NetCDF files already present.")
} else {
  message("Submitting ", length(jobs), " monthly jobs (sequential).")
  for (k in seq_along(jobs)) {
    j <- jobs[[k]]
    message(sprintf("[%d/%d] %s", k, length(jobs), j$target))
    download_one(user = cds_user, request = j$req, path = out_nc_dir, tries = 6, base_sleep = 20)
    
    fp <- file.path(out_nc_dir, j$target)
    if (!file.exists(fp) || !verify_nc(fp)) {
      stop("Downloaded file failed verification: ", fp)
    }
    
    # polite pause to reduce queue pressure
    Sys.sleep(8)
  }
}

# ----------------------------
# Process to daily means (per site) — robust version
# ----------------------------
process_site <- function(site,
                         drop_leap = FALSE,
                         allowed_ext = c("nc", "grib", "grb")) {
  
  # Match monthly files like: Oensingen_2014_01.nc / .grib / .grb
  ext_re <- paste(allowed_ext, collapse = "|")
  pat <- paste0("^", site, "_\\d{4}_\\d{2}\\.(", ext_re, ")$")
  files <- list.files(out_nc_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  # Sort by YYYYMM from filename
  key <- as.integer(gsub("^.*_(\\d{4})_(\\d{2})\\..*$", "\\1\\2", basename(files)))
  files <- files[order(key)]
  
  out <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    f  <- files[i]
    bn <- basename(f)
    yy <- as.integer(sub("^.*_(\\d{4})_(\\d{2})\\..*$", "\\1", bn))
    mm <- as.integer(sub("^.*_(\\d{4})_(\\d{2})\\..*$", "\\2", bn))
    
    # Open monthly stack of hourly total_cloud_cover
    r  <- suppressWarnings(terra::rast(f))
    nl <- terra::nlyr(r)
    if (nl == 0) stop("Empty raster: ", f)
    
    # Expected hours in month
    d1 <- as.Date(sprintf("%04d-%02d-01", yy, mm))
    d2 <- lubridate::ceiling_date(d1, "month") - lubridate::days(1)
    nd <- as.integer(lubridate::day(d2))
    expected <- nd * 24L
    
    # If CDS delivered fewer/more layers, trim to full 24-hour days
    if (nl < 24L) stop("Too few layers (<24) in ", f, " (have ", nl, ").")
    if (nl %% 24L != 0L) {
      nd_trim <- floor(nl / 24L)
      warning(sprintf("Layer count mismatch for %s: have %d, expected %d. Using first %d days.",
                      bn, nl, expected, nd_trim))
      r  <- r[[1:(nd_trim * 24L)]]
      nd <- nd_trim
      d2 <- d1 + (nd - 1L)
    }
    
    # Group 24-hour chunks per day and average
    idx   <- rep(seq_len(nd), each = 24L)
    r_day <- terra::tapp(r, index = idx, fun = "mean")
    
    vals  <- terra::values(r_day, dataframe = TRUE)
    dates <- seq(d1, d2, by = "1 day")
    ccov  <- as.numeric(vals[[1]])  # single variable
    
    out[[i]] <- tibble::tibble(
      sitename = site,
      date     = dates,
      ccov     = ccov
    )
    
    # Free memory early
    rm(r, r_day, vals); gc()
  }
  
  df <- dplyr::bind_rows(out) |> dplyr::arrange(date)
  if (drop_leap) {
    df <- df |>
      dplyr::filter(!(lubridate::month(date) == 2 & lubridate::mday(date) == 29))
  }
  df
}

# Run processing
dl <- lapply(sites$sitename, process_site)
dl <- dl[!vapply(dl, is.null, logical(1))]
if (length(dl) == 0) stop("No processed cloud cover data produced.")

# Save per-site and combined CSVs
invisible(lapply(dl, function(df_site) {
  out_site <- file.path(out_csv_dir, paste0(unique(df_site$sitename), "_ccov_daily.csv"))
  readr::write_csv(df_site, out_site)
  message("Saved: ", out_site, " (", nrow(df_site), " rows)")
}))

ccov_all <- dplyr::bind_rows(dl) |> dplyr::arrange(sitename, date)
out_all  <- file.path(out_csv_dir, "ccov_daily_all_sites.csv")
readr::write_csv(ccov_all, out_all)
message("Saved combined: ", out_all, " (", nrow(ccov_all), " rows)")