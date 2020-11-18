### lifetime simulation without treatment -------------------

### this can be done using the same functions as for internal validation
### only need to change the stopping expression & cleaned prep files

rm(list = ls())

dir_path <- "K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\"
source(paste0(dir_path, "src\\config\\define_dir_Rbase.R"))
source(paste0(dir_path, "src\\extrapolation\\preparation_functions.R"))
source(paste0(dir_path, "src\\model_paper\\extrapolation_functions_app_no_treat.R"))

### AKI
aki_without_antihyp <- c(0, 5.5, 17.9, 46.1, 0, 0, 0)/100000 
aki_with_antihyp <- aki_without_antihyp
aki_incidence <- rbind(aki_without_antihyp, aki_with_antihyp)
aki_qol_decr <- 0.0747

### stopping expression
stop_expr <- expression(floor(101 - (vX_i["age_cent"] * 10 + 64)))

### files
states_prefix <- "states"
mx_b_prefix <- "mx_b_"
mx_t_prefix <- "mx_t_" 
p_b_prefix <- "p_b_psa_app_no_treat_" 
mx_ckd_prefix <- "mx_ckd_"
alb_prefix <- "alb_"
covmeans_prefix <- "covmeans"
Q_prefix <- "Q"

### simulation parameters
output_file_prefix <- "df_output_no_treat"

### Estimation -------------------

n_cores <- detectCores() - 1
cf_prefix <- "2019-10-14_coeffs"
output_dir <- (paste0(dir_path, "output\\lifetime_datasets\\"))

### wrapper 
wrapper_est(
  data_for_extrapolation_dir = data_for_extrapolation_dir,
  data_risk_equations_dir = data_risk_equations_dir,
  output_dir = output_dir,
  mx_b_prefix = mx_b_prefix, mx_t_prefix = mx_t_prefix, p_b_prefix = p_b_prefix,
  state_prefix = state_prefix, 
  mx_ckd_prefix = mx_ckd_prefix, alb_prefix = alb_prefix, Q_prefix = Q_prefix,
  stop_expr = stop_expr,
  aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr, 
  output_file_prefix = output_file_prefix,
  cf_prefix = cf_prefix, n_cores = n_cores
) 

### PSA -------------------

n_cores <- 25
n_sim_start <- 701
n_sim_end <- 1000
boot_prefix <- "2019-09-05_boot_coeffs"
output_dir <- (paste0(dir_path, "output\\lifetime_datasets\\psa\\app\\no_treat\\"))

### wrapper
wrapper_psa(
  data_for_extrapolation_dir = data_for_extrapolation_dir,
  data_risk_equations_dir = data_risk_equations_dir,
  output_dir = output_dir,
  mx_b_prefix = mx_b_prefix, mx_t_prefix = mx_t_prefix, p_b_prefix = p_b_prefix,
  state_prefix = state_prefix,
  mx_ckd_prefix = mx_ckd_prefix, alb_prefix = alb_prefix, Q_prefix = Q_prefix,
  stop_expr = stop_expr,
  aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr,
  output_file_prefix = output_file_prefix,
  boot_prefix = boot_prefix, n_sim_start = n_sim_start, n_sim_end = n_sim_end, n_cores = n_cores
)