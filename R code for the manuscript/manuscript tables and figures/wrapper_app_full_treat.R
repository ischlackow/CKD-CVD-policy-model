### lifetime simulation with full treatment -------------------

rm(list = ls())

dir_path <- "K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\"
source(paste0(dir_path, "src\\config\\define_dir_Rbase.R"))
source(paste0(dir_path, "src\\extrapolation\\preparation_functions.R"))
source(paste0(dir_path, "src\\model_paper\\extrapolation_functions_app_full_treat.R"))

### stopping expression

stop_expr <- expression(floor(101 - (vX_i["age_cent"] * 10 + 64)))

### files

states_prefix <- "states"
mx_b_prefix <- "mx_b_"
mx_t_prefix <- "mx_t_" 
p_b_prefix <- "p_b_psa_app_no_treat_"  # is still suitable for the full treatment
mx_ckd_prefix <- "mx_ckd_"
alb_prefix <- "alb_"
covmeans_prefix <- "covmeans"
Q_prefix <- "Q"

### further files specific for treatment uptake scenarios

df_tx_flag_prefix  <- "df_tx_flag_antihyp_10_"
ckd_clin_prefix <- "ckd_clin_"
hyp_bp_prefix <- "hyp_bp_"
update_diab_prefix <- "update_diab_"

### parameters specific for treatment update scenarios

statin_diab <- c(0.001, 0.003)
statin_nvd <- 1.6e-06
statin_qol_decr <- 2.2438e-06

### AKI

aki_without_antihyp <- c(0, 5.5, 17.9, 46.1, 0, 0, 0)/100000 
aki_with_antihyp <- c(0, 6.4, 21.1, 54.4, 0, 0, 0)/100000
aki_incidence <- rbind(aki_without_antihyp, aki_with_antihyp)
aki_qol_decr <- 0.0747

### simulation parameters
output_file_prefix <- "df_output_full_treat"

### Estimation -------------------

n_cores <- detectCores() - 1
cf_prefix <- "2019-10-14_coeffs"
output_dir <- (paste0(dir_path, "output\\lifetime_datasets\\"))

# 1st statin: A20
tx_statin_vd_a20 <- c(0.784827031,	0.880787926,	0.899721945,	0.809731124,	0.816546135,	1, 0.899721945)
tx_statin_stroke_a20 <- c(0.771873377,	0.753051185,	0.887346328,	0.798077083,	0.805257459,	1.105363984, 0.887346328)
tx_statin_mi_a20 <- c(0.645751236,	0.642651088,	0.74175898,	0.821414406,	0.827856661,	0.850542403, 0.74175898)

# 2nd statin: A80
tx_statin_vd_a80 <- c(0.745797805, 0.857561591,	0.879927531,	0.774540562,	0.782438757,	1, 0.879927531)
tx_statin_stroke_a80 <- c(0.730922866,	0.709402657,	0.865297402,	0.761066667,	0.769363444,	1.128923003, 0.865297402) 
tx_statin_mi_a80 <- c(0.588950753,	0.585529773,	0.696545869,	0.788089306,	0.795577601,	0.822044108, 0.696545869)

# ramipril 10
tx_antihyp_vd <- c(0.87, rep(0.80, 6))
tx_antihyp_stroke <- c(0.85, rep(0.81, 6))
tx_antihyp_mi <- c(0.79, rep(0.81, 6))
tx_antihyp_hf <- c(0.90, rep(0.75, 6))

# aspirin
tx_antip_vd <- rep(0.87, 7)
tx_antip_stroke <- rep(1.00, 7)
tx_antip_mi <- rep(0.68, 7)

tx <- list(
  tx_statin_vd_a20 = tx_statin_vd_a20, tx_statin_stroke_a20 = tx_statin_stroke_a20, tx_statin_mi_a20 = tx_statin_mi_a20,
  tx_statin_vd_a80 = tx_statin_vd_a80, tx_statin_stroke_a80 = tx_statin_stroke_a80, tx_statin_mi_a80 = tx_statin_mi_a80,
  tx_antihyp_vd = tx_antihyp_vd, tx_antihyp_stroke = tx_antihyp_stroke, tx_antihyp_mi = tx_antihyp_mi, tx_antihyp_hf = tx_antihyp_hf,
  tx_antip_vd = tx_antip_vd, tx_antip_stroke = tx_antip_stroke, tx_antip_mi = tx_antip_mi
)

wrapper_est(
  data_for_extrapolation_dir = data_for_extrapolation_dir,
  data_risk_equations_dir = data_risk_equations_dir,
  output_dir = output_dir,
  mx_b_prefix = mx_b_prefix, mx_t_prefix = mx_t_prefix, p_b_prefix = p_b_prefix,
  state_prefix = state_prefix, 
  mx_ckd_prefix = mx_ckd_prefix, alb_prefix = alb_prefix, Q_prefix = Q_prefix,
  df_tx_flag_prefix  = df_tx_flag_prefix, ckd_clin_prefix = ckd_clin_prefix, hyp_bp_prefix = hyp_bp_prefix,
  update_diab_prefix = update_diab_prefix,
  stop_expr = stop_expr,
  statin_diab = statin_diab, statin_nvd = statin_nvd, statin_qol_decr = statin_qol_decr,
  aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr, 
  output_file_prefix = output_file_prefix,
  cf_prefix = cf_prefix, tx = tx, n_cores = n_cores
)

### PSA -----------------------------

n_cores <- 30
n_sim_start <- 1
n_sim_end <- 700
boot_prefix <- "2019-09-05_boot_coeffs"
output_dir <- (paste0(dir_path, "output\\lifetime_datasets\\psa\\app\\full_treat\\"))


wrapper_psa(
  data_for_extrapolation_dir = data_for_extrapolation_dir,
  data_risk_equations_dir = data_risk_equations_dir,
  output_dir = output_dir,
  mx_b_prefix = mx_b_prefix, mx_t_prefix = mx_t_prefix, p_b_prefix = p_b_prefix,
  state_prefix = state_prefix, 
  mx_ckd_prefix = mx_ckd_prefix, alb_prefix = alb_prefix, Q_prefix = Q_prefix,
  df_tx_flag_prefix  = df_tx_flag_prefix, ckd_clin_prefix = ckd_clin_prefix, hyp_bp_prefix = hyp_bp_prefix,
  update_diab_prefix = update_diab_prefix,
  stop_expr = stop_expr,
  statin_diab = statin_diab, statin_nvd = statin_nvd, statin_qol_decr = statin_qol_decr,
  aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr, 
  output_file_prefix = output_file_prefix,
  boot_prefix = boot_prefix, n_sim_start = n_sim_start, n_sim_end = n_sim_end, n_cores = n_cores
) 