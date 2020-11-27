###################################################################################
### internal validation
###################################################################################

rm(list = ls())

dir_path <- "K:\\SHARP_EE\\Monitoring CKD CHF\\master\\"
source(paste0(dir_path, "src\\config\\define_dir_Rbase.R"))
source(paste0(dir_path, "src\\extrapolation\\preparation_functions.R"))
source(paste0(dir_path, "src\\model_paper\\extrapolation_functions_intval.R"))

### stopping expression
stop_expr <- expression(6)

### files
states_prefix <- "states"
mx_b_prefix <- "mx_b_int_val_"
mx_t_prefix <- "mx_t_int_val_" 
p_b_prefix <- "p_b_int_val_psa_" 
mx_ckd_prefix <- "mx_ckd_int_val_"
alb_prefix <- "alb_int_val_"
covmeans_prefix <- "covmeans"
Q_prefix <- "Q"

### simulation parameters
output_file_prefix <- "df_output_int_val"

### Estimation -------------------

n_cores <- detectCores() - 1
cf_prefix <- "2019-10-14_coeffs"
output_dir <- (paste0(dir_path, "output\\lifetime_datasets\\"))

wrapper_est(
  data_for_extrapolation_dir = data_for_extrapolation_dir,
  data_risk_equations_dir = data_risk_equations_dir,
  output_dir = output_dir,
  mx_b_prefix = mx_b_prefix, mx_t_prefix = mx_t_prefix, p_b_prefix = p_b_prefix,
  state_prefix = state_prefix, 
  mx_ckd_prefix = mx_ckd_prefix, alb_prefix = alb_prefix, Q_prefix = Q_prefix,
  stop_expr = stop_expr,
  output_file_prefix = output_file_prefix,
  cf_prefix = cf_prefix, n_cores = n_cores
) 

### PSA -----------------------------
# 
# n_cores <- detectCores() - 1
# n_sim_start <- 1
# n_sim_end <- 1000
# boot_prefix <- "2019-09-05_boot_coeffs"
# output_dir <- (paste0(dir_path, "output\\lifetime_datasets\\intval\\"))
#  
# wrapper_est(
#    data_for_extrapolation_dir = data_for_extrapolation_dir,
#    data_risk_equations_dir = data_risk_equations_dir,
#    output_dir = output_dir,
#    mx_b_prefix = mx_b_prefix, mx_t_prefix = mx_t_prefix, p_b_prefix = p_b_prefix,
#    state_prefix = state_prefix, 
#    mx_ckd_prefix = mx_ckd_prefix, alb_prefix = alb_prefix, Q_prefix = Q_prefix,
#    stop_expr = stop_expr,
#    output_file_prefix = output_file_prefix,
#    boot_prefix = boot_prefix, n_sim_start = n_sim_start, n_sim_end = n_sim_end, n_cores = n_cores
#    ) 