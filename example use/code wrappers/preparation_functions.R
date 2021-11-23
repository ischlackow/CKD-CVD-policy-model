### states operations ##################################

# possible starting CV states
# 1 (no NFMVE) at the end of first CKD model run
.get_Ns_CV <- function(cycle) {
  if (cycle == 1)
    Ns <- 1
  if (cycle == 2)
    Ns <- c(1, 2, 3)
  return(Ns)
}

.get_Ns_CKD <- function(cycle) {
  if (cycle %in% 1:4)
    Ns <- 1:5
  if (cycle == 5)
    Ns <- 1:7
  if (cycle %in% 6:7)
    Ns <- 6:7
  return(Ns)
}

.get_next_N_CV <- function(N_CV){
  if (N_CV==1)
    Ns_CV1 <- c(1,2,3) else
      if(N_CV==2) 
        Ns_CV1 <- c(3,4) else
          if(N_CV==3)
            Ns_CV1 <- c(5) else
              if(N_CV==4)
                Ns_CV1 <- c(3,4) else
                  Ns_CV1 <- c(5)
                return(Ns_CV1)
}

.get_next_N_CKD <- function(N_CKD) {
  if(N_CKD %in% c(1:4))
    Ns_CKD1 <- c(1:5) else
      if(N_CKD==5)
        Ns_CKD1 <- c(1:7) else
          Ns_CKD1 <- c(6:7)
        return(Ns_CKD1)
}

.get_stateN <- function(states_info, N_CV, N_CKD) {
  N <- which(sapply(states_info, "[", "N_CV") == N_CV 
               & sapply(states_info, "[", "N_CKD") == N_CKD)
  names(N) <- NULL
  return(N)
}

get_state_combined <- function(N_CV, N_CKD) {
  return(paste(N_CV, N_CKD, sep = "_"))
}

### generate initial matrices ##################################

.clean_mx <- function(df, med_flag = FALSE) {
  
  # re-arrange
  df <- df %>%
    arrange(patid_temp) %>%
    select(patid_temp, everything()) %>%
    # remove redundant columns
    select(-c(patid, imd, bmi, bmi_cat, sbp_cent, chol_hdl)) %>%
    # re-define characteristics to align with the coefficients
    rename(bmi = bmi_cat_imputed, smoke = smoking_imputed, imd = imd_imputed,
           sbp_cent = sbp_cent_imputed, chol_ratio = chol_hdl_imputed) %>%
    # add a few more variables
    mutate(ckd_base = as.numeric(ckd)) %>%
    mutate(ckd3a = as.numeric(ckd == "ckd3a"), ckd3b = as.numeric(ckd == "ckd3b"), 
           ckd4 = as.numeric(ckd == "ckd4"), ckd5 = as.numeric(ckd == "ckd5")) %>%
    select(-ckd) %>%
    mutate(diab1 = as.numeric(diab_84 == 1), diab2 = as.numeric(diab_84 == 2)) %>%
    select(-diab_84) %>%
    mutate(baseline_CV = secondary, # todo: this could be done outside
           ckd_clin = as.numeric(ckd_base > 1 | acr > 1 | crd == 1)) %>% # todo: and this
    select(-c(smoking))
  # baseline medications
  if (med_flag == FALSE) 
    df <- mutate(df, statin_84 = 0, antihyp = 0, apc_base_84 = 0)
  
  return(df)
  
}

.gen_mm_b <- function(df, vars, vars_to_fact) {
  for (var in vars_to_fact){
    lab <- var[["lab"]]
    df[, lab] <- factor(df[, lab])
    if (var[["default"]] == FALSE)
      df[, lab] <- relevel(df[, lab], ref = var[["ref"]])
  }
  form <- as.formula(paste("~", paste(vars, collapse = "+"), sep = ""))
  mx_b <- data.frame(model.matrix(form, data = df))
  
  # add other variables
  mx_b <- mx_b %>%
    rename("(Intercept)" = X.Intercept.) %>%
    mutate(patid_temp = as.numeric(as.character(df$patid_temp))) %>%
    select(patid_temp, everything())
  
  return(mx_b)
}

.gen_mm_t <- function(df, vars){
  
  mx_t <- df %>% 
    select(vars) %>%
    mutate(patid_temp = as.numeric(as.character(df$patid_temp)))
  
  return(mx_t)
}

.clean_df_p_b <- function(mx_b,
                          cols_endpts, cols_hf, cols_qaly, cols_cost, cols_alive) {
  p_b <- mx_b %>%
    select(c(patid_temp, ckd_base)) %>%
    mutate(cycle = 0)
  
  # CV endpoints
  for (var in cols_endpts) {
    p_b <- mutate(p_b, !!var := 0)
  }
  p_b <- mutate(p_b, endpt_noNFMVE = 1)
  # states
  labs <- unique(sapply(states_info, "[[", 1))
  for (var in labs) {
    p_b <- mutate(p_b, !!var := 0)
  }
  for (N_CKD in (1:length(states_CKD))) {
    ids <- which(p_b$ckd_base == N_CKD)
    if (length(ids) > 0)
      p_b[ids, paste(1, N_CKD, sep = "_")] <- 1
  }
  # remove redundant columns
  p_b <- select(p_b, -ckd_base)
  # other columns
  for (var in c(cols_hf, cols_qaly, cols_cost)) {
    p_b <- mutate(p_b, !!var := 0)
  }
  # proportion surviving
  for (var in cols_alive) {
    p_b <- mutate(p_b, !!var := 1)
  }
  return(p_b)
}

### misc ##################################

.split_df_4 <- function(df, df_patients, prefix, remove_sex = TRUE, remove_secondary = TRUE, remove_patid_temp = FALSE,
                        output_dir = data_for_extrapolation_dir, to_vector = FALSE, to_matrix = TRUE){
  df <- merge(df, df_patients)
  for (sex_t in c("M", "F"))
    for (secondary_t in c(0, 1)){
      df_t <- df %>% 
        filter(sex == sex_t & secondary == secondary_t) %>%
        arrange(patid_temp)
      if (remove_sex)
        df_t <- select(df_t, -sex)
      if (remove_secondary)
        df_t <- select(df_t, -secondary)
      if (remove_patid_temp)
        df_t <- select(df_t, -patid_temp)
      df_t <- select(df_t, -patid)
      if (to_matrix)
        df_t <- as.matrix(df_t)
      if (to_vector)
        df_t <- unlist(df_t, use.names = FALSE)
      name <- paste(prefix, sex_t, secondary_t, sep = "_")
      assign(name, df_t)
      setwd(output_dir)
      save(list = name, file = paste0(name, ".Rdata"))
    }
}