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
  library(rsofun)
})
# read rds forcing data file previously uploaded for cn model as drivers data
drivers <- readRDS("./data/FLX_CH-Oe2_FLUXNET2015_FULLSET_2004-2023_1-3/rsofun_driver_data_v3.4.2.rds") %>% 
  dplyr::filter(sitename == "CH-Oe2") 

# replace forcing data forcing_df <- drivers$forcing[[1]] with a new csv file

# read csv file with forcing data
forcing_df <- read_csv("./data/01_data_prep/09_forcing_data_2004-2023.csv") 

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
n_input <- read.csv("./data/mangement_data/04_23choe2_ninput.csv")

#fharv 
fharv <- read.csv("./data/mangement_data/04_23ch0e2_fharv.csv")
# convert to fraction from percentage to be consistent with the rest of the data
fharv$fharv <- fharv$fharv / 100

#seed 
seed <- read.csv("./data/mangement_data/04_23ch0e2_seed.csv")

#yieldcn  
yieldcn <- read.csv("./data/mangement_data/04_23ch0e2_yieldcn.csv")

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

# Define model parameters taken from p model script and other scripts
pars <- list(
  # Paramteres taken from P-model
  kphio = 0.045,             # setup ORG in Stocker et al. 2020 GMD  #typical ranges: 0.04-0.1 #lma / fapar / gpp / leaf c / cturnover / Vcmax 
  kphio_par_a = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b = 1.0,
  soilm_thetastar = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous # cturnover / Vcmax
  tau_acclim = 30.0,
  kc_jmax = 0.41,
  
  # Plant
  f_nretain = 0.500000, #0.5 # cturnover
  fpc_tree_max = 0.8000000, #typical range: 0.8-0.95. #fapar / c leaf
  growtheff = 0.6000000, #0.45  #typical range: 0.3-0.7 #lma / c leaf / cturnover
  r_root = 2*0.913000,
  r_sapw = 2*0.044000,
  exurate = 0.003000,
  
  k_decay_leaf = 1.90000, #1.90000, #typical range: 0.5 - 2.0 # cleaf / cturnover
  k_decay_root = 1.90000, # cturnover
  k_decay_labl = 1.90000, #1.90000, #typical range: 0.5 - 2.0 # cleaf / cturnover
  k_decay_sapw = 1.90000, # cturnover
  
  r_cton_root = 37.0000, # cturnover
  r_cton_wood = 100.000, # cturnover
  r_cton_seed = 15.0000, # cturnover
  nv_vcmax25 = 5000, #0.02 * 13681.77, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.695; previously: 5000.0, #lma / Vcmax
  ncw_min = 0.056, #0.08 * 1.116222, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.691; previously used: 0.056, #lma
  r_n_cw_v = 0.05, #0.1, # assumed that LMA is independent of Vcmax25; previously: 0.1, #lma /Vcmax #typical range: 0.05-0.2 
  r_ctostructn_leaf = 30.0000, #1.3 * 45.84125, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.699; previously used: 80.0000, #typical range: 40-80 #lma / fapar / cleaf 
  kbeer = 0.350000, #0.400000 #fapar #typical range: 0.3-0.8
  
  # Phenology (should be PFT-specific)
  gddbase = 5.0, #typical range: 5-10 #c leaf 
  ramp = 0.0,
  phentype = 2.0,
  
  # Soil physics (should be derived from params_soil, fsand, fclay, forg, fgravel)
  perc_k1 = 5.0,
  thdiff_wp = 0.2,
  thdiff_whc15 = 0.8,
  thdiff_fc = 0.4,
  forg = 0.01,
  wbwp = 0.029,
  por = 0.421,
  fsand = 0.82,
  fclay = 0.06,
  fsilt = 0.12,
  
  # Water and energy balance
  kA = 107,
  kalb_sw = 0.17,
  kalb_vis = 0.03,
  kb = 0.20,
  kc = 0.25,
  kCw = 1.05,
  kd = 0.50,
  ke = 0.0167,
  keps = 23.44,
  kWm = 220.0,
  kw = 0.26,
  komega = 283.0,
  maxmeltrate = 3.0,
  
  # Soil BGC
  klitt_af10 = 1.2,
  klitt_as10 = 0.35,
  klitt_bg10 = 0.35,
  kexu10 = 50.0,
  ksoil_fs10 = 0.021,
  ksoil_sl10 = 7.0e-04,
  ntoc_crit1 = 0.45,
  ntoc_crit2 = 0.76,
  cton_microb = 10.0,
  cton_soil = 9.77,
  fastfrac = 0.985,
  
  # N uptake
  eff_nup = 0.0001000,
  minimumcostfix = 1.000000,
  fixoptimum = 25.15000,
  a_param_fix = -3.62000,
  b_param_fix = 0.270000,
  
  # Inorganic N transformations (re-interpreted for simple ntransform model)
  maxnitr =  0.00005,
  
  # Inorganic N transformations for full ntransform model (not used in simple model)
  non = 0.01,
  n2on = 0.0005,
  kn = 83.0,
  kdoc = 17.0,
  docmax = 1.0,
  dnitr2n2o = 0.01,
  
  # Additional parameters - previously forgotten
  frac_leaf = 0.4,         # after wood allocation  #typical range: 0.4-0.8 #lma / fapar / c leaf / cturnover
  frac_wood = 0,           # highest priority in allocation  #typical range: 0-0.6 #lma / cleaf 
  frac_avl_labl = 0.1, #0.1   #typical range: 0.1-0.3 # c leaf 
  
  # for development
  tmppar = 9999,
  
  # simple N uptake module parameters
  nuptake_kc = 600,
  nuptake_kv = 5,
  nuptake_vmax = 0.3 #0.3 #typical range: 0.1-0.3 # cleaf / cturnover
)

# make dir
if (!dir.exists("output")) {
  dir.create("output")
}
# Function to run the model and save the output
cnmodel_run_save <- function(drivers, pars, save_path, file_basename) {
  output_04_23 <- runread_cnmodel_f(drivers = drivers, 
                                    par = pars,
                                    ncores = 12,
                                    makecheck = TRUE,
                                    parallel = TRUE
  )
  rds_filename <- paste0(file_basename, ".rds")
  if (!dir.exists(save_path)) {
    stop("Directory does not exist: ", save_path)
  }
  saveRDS(output_04_23, file = file.path(save_path, rds_filename))
  message("Model output saved: ", file.path(save_path, rds_filename))
  invisible(output_04_23)
}
# Example usage:
# Run the model and save output
cnmodel_run_save(drivers, pars, "output", "output_04_23")
# Read the output
output_04_23 <- readRDS("output/output_04_23.rds")
# FLUXNET data for site CH-Oe2 (Oensingen, Switzerland)
drivers$params_siml <- modify_params(drivers$params_siml, 15, 2)
drivers$params_siml[[1]]$c_only <- TRUE

# # make drivers$forcing[[1]]$netrad whole column NA
# drivers$forcing[[1]]$netrad <- NA

# save the updated drivers data
saveRDS(drivers, "./data/CH-Oe2_2004-2023_final_ready_for_CNmodel_run.rds")