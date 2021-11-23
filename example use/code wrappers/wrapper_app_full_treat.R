###############################################################
### lifetime simulation with full treatment -------------------
###############################################################

### define directories ----

rm(list = ls())

# PLEASE CHANGE THESE LINES TO YOUR DIRECTORY
project_dir <- "Z:/HERC/CKD-CVD+parameter+data/CKD-CVD-policy-model/example use"
data_for_extrapolation_dir <- file.path(project_dir, "input patient data")
data_risk_equations_dir <- file.path(project_dir, "risk equations")

### load functions ----

# misc functions
source(file.path(project_dir, "code wrappers/preparation_functions.R"))
# functions for extrapolation
source(file.path(project_dir, "code wrappers/extrapolation_functions_app_full_treat.R"))

### model parameters ----

### stopping expression ###

# run until each patient has reached 101 years old
stop_expr <- expression(floor(101 - (vX_i["age_cent"] * 10 + 64)))
# another example: run for 5 years
# stop_expr <- expression(6) # need 6 here as start counting from 0, ie baseline

### filename prefixes ###

# Note that each patient-related file prefix will be associated with 4 files,
# ie one for each gender (F/M) & type of secondary prevention (0/1) combination

### files for CV risk equations

# matrix of baseline characteristics (for CV risk equations)
mx_b_prefix <- "mx_b_"
# matrix with time-updated characteristics (for CV risk equations)
mx_t_prefix <- "mx_t_" 

### files for CKD risk equations

# matrix with values for the CKD risk equations
mx_ckd_prefix <- "mx_ckd_"
# data covariate means for CKD risk equations
covmeans_prefix <- "covmeans"
# allowed transitions in the CKD model
Q_prefix <- "Q"

### general model files

# description  / parameters of states and endpoints
states_prefix <- "states"
# starting (ie baseline) value of endpoints and state distribution
p_b_prefix <- "p_b_ce_"

# files specific for treatment uptake scenarios

# matrix with treatment values for statins/antihypertensives depending on cycle (ie patient's age)
df_tx_flag_prefix  <- "df_tx_flag_"
# matrix with patients' baseline clinical CKD status
ckd_clin_prefix <- "ckd_clin_"
# matrix with patients' baseline hypertension status
hyp_bp_prefix <- "hyp_bp_"
# matrix with patients' baseline ACR status
alb_prefix <- "alb_"
# matrix indicating whether patient's diabetes status should be updated (eg because they don't have baseline diabetes)
update_diab_prefix <- "update_diab_"

### adverse events

# annual probability of developing statin-induced diabetes
# values for statin 1 (A20) & statin 2 (A80)
statin_diab <- c(0.00098, 0.00299)
# annual probability of statin-induced nonvascular mortality
statin_nvd <- 0.0000016
# QoL decrement associated with taking statins due to developing myopathy 
statin_qol_decr_myopathy <- 1.53699e-07
# QoL relative decrement associated with taking statins due to developing rhabdomyolysis
statin_qol_scale_rhabdo <- 1 * (1 - 3.4/100000) + 0.0267 * (3.4 / 100000)
# cost decrement associated with statin adverse effects
statin_cost_decr <- 0.1344025

### acute kidney injury (AKI)

# annual probability of developing AKI if not on anti-hypertensive drugs (ie due to CKD)
# for stages 2, 3a, 3b, 4, 5, transplant, dialysis
aki_without_antihyp <- c(0, 5.5, 17.9, 46.1, 0, 0, 0)/100000 
# annual probability of developing AKI if on anti-hypertensive drugs (ie due to CKD and anti-hypertensives)
# for stages 2, 3a, 3b, 4, 5, dialysis, transplant
aki_with_antihyp <- c(0, 6.4, 21.1, 54.4, 0, 0, 0)/100000
aki_incidence <- rbind(aki_without_antihyp, aki_with_antihyp)
# QoL decrement associated with AKI 
aki_qol_decr <- 0.0747
# cost decrement associated with AKI 
aki_cost_decr <- 1994.56

### treatment costs (June 2019)

# statins
t1_price <- (0.81 / 28) * 365.25
t2_price <- (1.65 / 28) * 365.25
statin_price <- c(t1_price, t2_price)

# antihypertensives
ram10_price <- (1.08 / 28) * 365.25
antihyp_price <- ram10_price

# antiplatelets
asp75_price <- (0.71 / 28) * 365.25
antip_price <- asp75_price

### discounting

disc <- 0.035

### output filename

output_file_prefix <- "df_output_full_treat"

### Estimation -------------------

# number of cores to run the model
n_cores <- detectCores() - 1

# coefficient filename
cf_prefix <- "2019-10-14_coeffs"

# output directory
output_dir <- file.path(project_dir, "output/estimate")

### treatment effects ###

# effects are presented by CKD stage in order ckd 2; ckd 3a; ckd 3b; ckd 4; ckd 5; transplant; dialysis

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

### run the wrapper ###

source(file.path(project_dir, "code wrappers/extrapolation_functions_app_full_treat.R"))

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
  statin_diab = statin_diab, statin_nvd = statin_nvd, 
  statin_qol_decr_myopathy = statin_qol_decr_myopathy,
  statin_qol_scale_rhabdo = statin_qol_scale_rhabdo, statin_cost_decr = statin_cost_decr,
  aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr, aki_cost_decr = aki_cost_decr,
  statin_price = statin_price, antihyp_price = antihyp_price, antip_price = antip_price,
  output_file_prefix = output_file_prefix,
  cf_prefix = cf_prefix, tx = tx, n_cores = n_cores, disc = disc
)

### PSA -----------------------------

# number of cores to run the model
n_cores <- detectCores() - 1

# simulation start & end number
n_sim_start <- 1
n_sim_end <- 5

# coefficient filename
boot_prefix <- "2019-09-05_boot_coeffs"

# output directory
output_dir <- file.path(project_dir, "output/psa")

### run the wrapper ###

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
  statin_diab = statin_diab, statin_nvd = statin_nvd, 
  statin_qol_decr_myopathy = statin_qol_decr_myopathy,
  statin_qol_scale_rhabdo = statin_qol_scale_rhabdo, statin_cost_decr = statin_cost_decr,
  aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr, aki_cost_decr = aki_cost_decr,
  statin_price = statin_price, antihyp_price = antihyp_price, antip_price = antip_price,
  output_file_prefix = output_file_prefix,
  boot_prefix = boot_prefix, n_cores = n_cores, disc = disc,
  n_sim_start = n_sim_start, n_sim_end = n_sim_end
)