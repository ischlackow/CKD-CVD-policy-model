###################################################################
#                       libraries              	                  #
###################################################################

library(tidyverse)
library(parallel)
library(doParallel)

###################################################################
#                        preliminary calculations        		      #
###################################################################

### risk equations

prelimRiskEquations <- function(vars_b, mx_b, coeff_NVD, coeff_VD, 
                                coeff_VD_Stroke, coeff_MI_Stroke_VD, coeff_HF) {
  
  ### split coefficients into baseline and time-updated
  for (var in c("coeff_NVD", "coeff_VD", "coeff_VD_Stroke", "coeff_MI_Stroke_VD", "coeff_HF")) {
    names_b <- intersect(vars_b, names(get(var)))
    names_t <- setdiff(names(get(var)), names_b)
    assign(paste(var), append(get(var), 
                              list(coeff_b = get(var)[names_b],
                                   coeff_t = get(var)[names_t])))
  }
  
  ### produce baseline lambdas
  lam_NVD_b_All <- exp(.lamM(beta = coeff_NVD$coeff_b, X = mx_b))
  lam_VD_b_All <- exp(.lamM(beta = coeff_VD$coeff_b, X = mx_b))
  lam_VD_Stroke_b_All <- exp(.lamM(beta = coeff_VD_Stroke$coeff_b, X = mx_b))
  lam_MI_Stroke_VD_b_All <- exp(.lamM(beta = coeff_MI_Stroke_VD$coeff_b, X = mx_b))
  lam_HF_b_All <- exp(.lamM(beta = coeff_HF$coeff_b, X = mx_b))
  
  return(list(
    lambdas = list(lam_NVD_b_All = lam_NVD_b_All, lam_VD_b_All = lam_VD_b_All, lam_VD_Stroke_b_All = lam_VD_Stroke_b_All,
                   lam_MI_Stroke_VD_b_All = lam_MI_Stroke_VD_b_All, lam_HF_b_All=lam_HF_b_All),
    coeffs = list(coeff_NVD = coeff_NVD, coeff_VD = coeff_VD, coeff_VD_Stroke = coeff_VD_Stroke, 
                  coeff_MI_Stroke_VD=coeff_MI_Stroke_VD, coeff_HF=coeff_HF)))
}

prelimCost <- function(vars_b, mx_b, coeff_cost){
  names_b <- intersect(vars_b, names(get("coeff_cost")))
  names_t <- setdiff(names(get("coeff_cost")[["coeff"]]), names_b)
  assign(paste("coeff_cost"), append(get("coeff_cost"), 
                                     list(coeff_b = get("coeff_cost")[names_b],
                                          coeff_t = get("coeff_cost")[names_t])))
  lam_cost_b <- .lamM(beta = coeff_cost$coeff_b, X = mx_b)
  return(list(lambdas = list(lam_cost_b=lam_cost_b),
              coeff = list(coeff_cost = coeff_cost)))
}

prelimQol <- function(vars_b, mx_b, coeff_qol){
  names_b <- intersect(vars_b, names(get("coeff_qol")))
  names_t <- setdiff(names(get("coeff_qol")), names_b)
  assign(paste("coeff_qol"), append(get("coeff_qol"), 
                                    list(coeff_b = get("coeff_qol")[names_b],
                                         coeff_t = get("coeff_qol")[names_t])))
  lam_qol_b <- .lamM(beta = coeff_qol$coeff_b, X = mx_b)
  return(list(lambdas = list(lam_qol_b=lam_qol_b),
              coeff = list(coeff_qol = coeff_qol)))
}

###################################################################
#                  extracting probabilities            		        #
###################################################################

# calculate \beta*X
.lamVec <- function(beta, X) {
  X <- X[names(beta)]
  return(X %*% beta)
}

# internal function for calculating beta * X
.lamM <- function(beta, X) {
  X <- X[, names(beta)] 
  #if (length(beta) == 1)
  #X <- as.matrix(X, ncol = length(beta))
  return(X %*% beta)
}

# cumulative hazard for Weibull distribution
H_w <- function(lam_b, beta, X, p, t) {
  lam <- lam_b * exp(.lamVec(beta, X)) 
  return(lam * ((t-1)^exp(p)-t^(exp(p))))
}

# cumulative hazard for Gompertz distribution
H_g <- function(lam_b, beta, X, g, t) {
  lam <- lam_b * exp(.lamVec(beta, X))
  return((lam / g) * exp(g * t) * (- 1 + exp(-g))) 
}


# Calculate the transition probability from a survival model 
tp_survival_w <- function(lam_b, X, beta, anc, t) {
  H1 <- H_w(lam_b = lam_b, beta = beta, X = X, p = anc, t = t)
  return(1 - exp(H1))
}

tp_survival_g <- function(lam_b, X, beta, anc, t) {
  H1 <- H_g(lam_b = lam_b, beta = beta, X = X, g = anc, t = t)
  return(1 - exp(H1))
}

get_pcond <- function(p, pNVD) {
  return(p * (1 - pNVD))
}

get_tp_mx_ckd <- function(v_ckd, t, covmeans, Q, q_mx, tp_esrd){
  
  # read off covariate values
  AGE60 <- v_ckd["AGE60"] + t
  covlist <- c(AGE60, v_ckd["SEX"], v_ckd["HF"], v_ckd["CA"])
  # Evaluate intensity matrix at the covariate values
  logest <- q_mx[["Q_logbaseline"]]
  for (j in 1:4) 
    logest <- logest + q_mx[[j + 1]] * (covlist[j] - covmeans[j])
  res <- exp(logest)
  res[Q == 0] <- 0
  diag(res) <- 0
  diag(res) <- -rowSums(res)
  # turn into a TP matrix. DO NOT APPLY MISSCLASSIFICATION MATRIX
  temp <- expm::expm(res)
  
  # re-distribute death across other states
  # but note only one diagonal below the main diagonal is permitted so keep it that way
  temp[1:3, 1:5] <- temp[1:3, 1:5] + temp[1:3, 6] / 5
  temp[4, 2:5] <- temp[4, 2:5] + temp[4, 6] / 4
  temp[5, 3:5] <- temp[5, 3:5] + temp[5, 6] / 3
  
  # add to the matrix with ESRD-related transitions
  age_index <- (AGE60 >= -10) + 1
  tp_mx <- tp_esrd[[age_index]]
  tp_mx[1:5, 1:5] <- temp[1:5, 1:5]
  tp_mx[5, 1:5] <- tp_mx[5, 1:5] * (1 - sum(tp_mx[5, 6:7]))
  
  # save & return
  return(tp_mx)
}

cycle_cost <- function(X, coeff_cost_t, j, lam_cost_b){
  time_costs <- .lamVec(beta = coeff_cost_t, X = X)
  costs <- lam_cost_b + time_costs
  return(drop(costs))
}

cycle_qol <- function(X, coeff_qol_t, j, lam_qol_b){
  # absolute decrements
  decr_states <- c("ckd3a", "ckd3b", "ckd4", "ckd5", "dialysis", "transplant", "age_cent", "diab1", "diab2")
  decr_qol <- .lamVec(beta = coeff_qol_t[decr_states], X = X)
  # multipliers
  event_states <- c("endpt_noNFMVE", "endpt_MI", "endpt_Stroke", "endpt_prevMI", "endpt_prevStroke")
  event_states_HF <- c("HF_hosp", "prev_HF")
  event_qol <- .lamVec(beta = coeff_qol_t[event_states], X = X) 
  event_qol_HF <- .lamVec(beta = coeff_qol_t[event_states_HF], X = X) + (1 - sum(X[event_states_HF]))
  qol <- (lam_qol_b + decr_qol) * event_qol * event_qol_HF
  return(drop(qol))
}

gen_nextCycle <- function(j, secondary, states_info, endpts_CV, names_mx_mon,
                          tp_survival_nvd, coeff_NVD_t, anc_nvd, lam_NVD_b, 
                          tp_survival_vd, coeff_VD_t, anc_vd, lam_VD_b,
                          tp_survival_vd_stroke, coeff_VD_Stroke_t, anc_vd_stroke, lam_VD_Stroke_b,
                          tp_survival_vd_stroke_mi, coeff_MI_Stroke_VD_t, anc_vd_stroke_mi, lam_MI_Stroke_VD_b,
                          tp_survival_hf, coeff_HF_t, anc_hf, lam_HF_b,
                          tx_statin, tx_antihyp, tx_antip,
                          coeff_qol_t, lam_qol_b, coeff_names_cost_qol,
                          vP_b, vX_i, 
                          covmeans_i, Q_ckd_i, mx_ckd_i, tp_esrd, Q,
                          treat_i, ckd_clin_i, hyp_bp_i, update_diab_i,
                          j_G2, j_G3a, j_G3b, j_G4, j_G5,
                          statin_diab, statin_nvd, statin_qol_decr,
                          aki_incidence, aki_qol_decr){
  
  ### update relevant covariates ----------------------------------
  
  # last year counter
  j_prev <- j - 1
  
  # read off treatment values for the patient in 1_1 at the start of the year
  treat_i_j_start <- treat_i[treat_i[, "cycle"] == j_prev, ]
  need_statin_i_j_start <- as.numeric(treat_i_j_start[ "need_statins"])
  need_antihyp_i_j_start <- as.numeric(treat_i_j_start["need_antihyp"])
  
  ### update output probability vectors
  
  # 1st output vector
  vP_1 <- vP_b
  vP_1["cycle"] <- j
  # carry over the count of HF hospitalisations
  # just add up HF_hosp & prev_HF and scale by those alive at the start of the year
  # proportion alive at the end of the last year
  n_alive_at_start <- vP_b["alive_at_end"]
  vP_1["alive_at_start"] <- n_alive_at_start
  #prev_HF <- as.numeric(vP_b["HF_hosp"] * vP_1["alive_last"] + vP_b["prev_HF"])
  d_last_year <- vP_1["endpt_VD"] + vP_1["endpt_NVD"]
  prev_HF <- as.numeric(vP_b["prev_HF"] * (n_alive_at_start / (n_alive_at_start + d_last_year)) + vP_b["HF_hosp"]) 
  vP_1["prev_HF"] <- prev_HF
  
  # 2nd output vector - after the CV model is run
  vP_2 <- vP_1
  
  # generate vector with time-updated covariates
  vX_i["age_cent"] <- vX_i["age_cent"] + 0.1 * j_prev
  
  ### run CKD model ----------------------------------
  
  tp_mx_ckd <- get_tp_mx_ckd(v_ckd = mx_ckd_i, t = j_prev, covmeans = covmeans_i, Q = Q, q_mx = Q_ckd_i, tp_esrd = tp_esrd)
  
  vP_1[j_G2] <- vP_b[j_G2] %*% tp_mx_ckd
  vP_1[j_G3a] <- vP_b[j_G3a] %*% tp_mx_ckd
  vP_1[j_G3b] <- vP_b[j_G3b] %*% tp_mx_ckd
  vP_1[j_G4] <- vP_b[j_G4] %*% tp_mx_ckd
  vP_1[j_G5] <- vP_b[j_G5] %*% tp_mx_ckd
  
  ### run CV model ----------------------------------
  
  ### extract possible end states and outcomes
  
  # states0 - only needed to extract states2 for cycle 1 (otherwise states2 are just all states)
  # states1 - after the CKD model has been run, ie change the CKD component
  # states2 - after the CV model has been run, ie change the CV component; only need the lab
  if (j == 1) {
    states1 <- states_info[unlist(1:7)]
    states1_lab <- names_mx_mon[1:7]
    l1 <- 7
    states0 <- states_info[unlist(1:5)]
    states2_num <- unique(unlist(lapply(states0, '[[', 'Ns2'), use.names = FALSE))
    states2 <- states_info[unlist(states2_num)]
    states2_lab <-  sapply(states2, '[[', 'lab')
  } else {
    states1 <- states_info
    states1_lab <- names_mx_mon
    l1 <- 35
    states2_lab <- names_mx_mon
  }
  
  mx_cv_rownames <- states1_lab
  mx_cv_colnames <- c(endpts_CV, states2_lab, c("HF_hosp", "qaly"))
  
  mx_cv <- matrix(0, nrow = l1, ncol = length(mx_cv_colnames),
                  dimnames = list(mx_cv_rownames, mx_cv_colnames))
  
  ### loop through each starting state  ----------------------------------
  
  #s <- 1
  #s <- 5
  for (s in 1:l1) {
    
    ### read off state characteristics  ####################
    
    state1 <- states1[[s]]
    lab <- state1[["lab"]]
    N_CV1 <- state1[["N_CV"]]
    N_CKD1 <- state1[["N_CKD"]]
    eligible_stroke <- state1[["eligible_stroke"]]
    eligible_MI <- state1[["eligible_MI"]]
    
    # covariates
    covars_ckd_5 <- state1[["covars_ckd_5"]]
    covars_ckd_all <- state1[["covars_ckd_all"]]
    covars_prev_event <- state1[["covars_prev_event"]]
    
    ### update vX_i 
    
    vX_i_s <- c(vX_i, covars_ckd_5, covars_prev_event)
    
    ### ckd indicator ####################
    
    ckd_state <- (N_CKD1 != 1)
    
    ### treatment indicators at start  ####################
    
    if (N_CV1 != 1 | secondary){
      # start in a secondary prevention state or come from secondary prevention population
      # statins: everyone is on A80
      on_statin_1_start <- 0
      on_statin_2_start <- 1
      # antihypertensives: everyone
      on_antihyp_start <- 1
      # antiplatelets: everyone
      on_antip_start <- 1
    } else {
      # start in a noNFMVE state & come from primary prevention population
      # statins: A20 for those with CKD starting state OR need statins at the start for whatever reason (takes baseline CKD into account)
      on_statin_1_start <- max(ckd_state, need_statin_i_j_start)
      on_statin_2_start <- 0
      # antihypertensives: everyone who needs statins for whatever reason OR ckd + stage 1 hypertension
      on_antihyp_start <- max(ckd_state * hyp_bp_i, need_antihyp_i_j_start)
      # antiplatelets: anyone
      on_antip_start <- 0
    }
    not_on_statin_start <- 1 - on_statin_1_start - on_statin_2_start
    not_on_antihyp_start <- 1 - on_antihyp_start
    not_on_antip_start <- 1 - on_antip_start
    
    ### extracting treatment effects ####################
    
    ### statins
    
    # weight the Tx effect
    tx_statin_s <- not_on_statin_start * rep(1, 3) + on_statin_1_start * tx_statin[["prim"]][, N_CKD1] + on_statin_2_start * tx_statin[["sec"]][, N_CKD1]
    
    # extract effect for each endpoint
    tx_statin_vd_s <- tx_statin_s[1]
    tx_statin_stroke_s <- tx_statin_s[2]
    tx_statin_mi_s <- tx_statin_s[3]
    
    # udpate diabetes
    if (update_diab_i)
      vX_i_s["diab2"] <- vX_i_s["diab2"] + on_statin_1_start * statin_diab[1] + on_statin_2_start * statin_diab[2]
    
    ### anti-hypertensives
    
    # weight the Tx effect
    tx_antihyp_s <- not_on_antihyp_start * rep(1, 4) + on_antihyp_start * tx_antihyp[, N_CKD1]
    
    # extract effect for each endpoint
    tx_antihyp_vd_s <- tx_antihyp_s[1]
    tx_antihyp_stroke_s <- tx_antihyp_s[2]
    tx_antihyp_mi_s <- tx_antihyp_s[3]
    tx_antihyp_hf_s <- tx_antihyp_s[4]
    
    ### antiplatelets
    
    # weight the Tx effect
    tx_antip_s <- not_on_antip_start * rep(1, 3) + on_antip_start * tx_antip[, N_CKD1]
    
    # extract effect for each endpoint
    tx_antip_vd_s <- tx_antip_s[1]
    tx_antip_stroke_s <- tx_antip_s[2]
    tx_antip_mi_s <- tx_antip_s[3]
    
    ### generate the starting vector ####################
    
    mx_cv_s <- mx_cv[s, ]
    
    ### state transitions ####################
    
    ### NVD
    
    pNVD <- do.call(tp_survival_nvd, args = list(lam_b = lam_NVD_b, X = vX_i_s, beta = coeff_NVD_t, anc = anc_nvd, t = j))
    # excess mortality from statins
    pNVD <- pNVD + (on_statin_1_start + on_statin_2_start) * statin_nvd
    
    # update the probability matrix
    mx_cv_s["endpt_NVD"] <- pNVD 
    
    ### VD
    
    pVD <- do.call(tp_survival_vd, args = list(lam_b = lam_VD_b, X = vX_i_s, beta = coeff_VD_t, anc = anc_vd, t = j))
    pVD_w <- pVD * tx_statin_vd_s * tx_antihyp_vd_s * tx_antip_vd_s
    pVD_w <- get_pcond(p = pVD_w, pNVD = pNVD)
    # update probability matrix
    mx_cv_s["endpt_VD"] <- pVD_w
    
    ### stroke
    
    if (eligible_stroke){
      pStroke_VD <- do.call(tp_survival_vd_stroke, args = list(lam_b = lam_VD_Stroke_b, X = vX_i_s, beta = coeff_VD_Stroke_t, anc = anc_vd_stroke, t = j))
      pStroke <- max(pStroke_VD - pVD, 0)
      pStroke <- pStroke * tx_statin_stroke_s * tx_antihyp_stroke_s * tx_antip_stroke_s 
      pStroke <- get_pcond(p = pStroke, pNVD = pNVD)
    } else pStroke <- 0
    
    mx_cv_s["endpt_Stroke"] <- pStroke
    mx_cv_s[get_state_combined(3, N_CKD1)] <- pStroke
    
    ### MI
    
    if(eligible_MI){
      pMI_Stroke_VD <- do.call(tp_survival_vd_stroke_mi, args = list(lam_b = lam_MI_Stroke_VD_b, X = vX_i_s, beta = coeff_MI_Stroke_VD_t,
                                                                     anc = anc_vd_stroke_mi, t = j))
      
      pMI <- if (pStroke_VD > pVD) max(pMI_Stroke_VD - pStroke_VD, 0) else
        max(pMI_Stroke_VD - pVD, 0)
      pMI <- pMI * tx_statin_mi_s * tx_antihyp_mi_s * tx_antip_mi_s 
      pMI <- get_pcond(p = pMI, pNVD = pNVD)
    } else pMI <- 0
    
    mx_cv_s["endpt_MI"] <- pMI
    mx_cv_s[get_state_combined(2, N_CKD1)] <- pMI
    
    ### remainder
    N_CV2 <- state1[["N_CV_remaining"]]
    p_d <- pNVD + pVD_w 
    p_remained <- 1 - sum(mx_cv_s[states2_lab]) - p_d
    mx_cv_s[get_state_combined(N_CV2, N_CKD1)] <- p_remained
    
    # previous event
    mx_cv_s["endpt_noNFMVE"] <- mx_cv_s[get_state_combined(1, N_CKD1)]
    if (j > 1) {
      mx_cv_s["endpt_prevMI"] <- mx_cv_s[get_state_combined(4, N_CKD1)]
      mx_cv_s["endpt_prevStroke"] <- mx_cv_s[get_state_combined(5, N_CKD1)]
    }
    
    ### HF 
    
    # update previous event variable
    vX_i_s["prev_event1"] <- sum(mx_cv_s[intersect(paste(c(2, 4), N_CKD1, sep = "_"), mx_cv_colnames)])
    vX_i_s["prev_event2"] <- sum(mx_cv_s[intersect(paste(c(3, 5), N_CKD1, sep = "_"), mx_cv_colnames)])
    
    pHF <- do.call(tp_survival_hf, args = list(lam_b = lam_HF_b, X = vX_i_s, beta = coeff_HF_t, anc = anc_hf, t = j))
    pHF <- pHF * tx_antihyp_hf_s
    mx_cv_s["HF_hosp"] <- (1 - prev_HF) * pHF
    
    ### probability of death
    
    d_midyear <- p_d / 2
    alive_midyear <- 1 - d_midyear
    
    ### qalys  ####################
    
    # vector with covariates
    X <- c(vX_i_s[c("age_cent", "diab1", "diab2")], covars_ckd_all, mx_cv_s[coeff_names_cost_qol], prev_HF = prev_HF)
    
    ### main effects
    
    qaly <- cycle_qol(X = X, coeff_qol_t = coeff_qol_t, j = j, lam_qol_b = lam_qol_b)
    
    ### adverse events decrements
    
    # due to rhabdo & myopathy
    qaly <- qaly - (on_statin_1_start + on_statin_2_start) * statin_qol_decr
    # due to AKI
    aki_risk <- (1 - on_antihyp_start) * aki_incidence[1, N_CKD1] + on_antihyp_start * aki_incidence[2, N_CKD1]
    qaly <- qaly - aki_risk * aki_qol_decr
    
    ### weight by 1/2 for those who died
    
    #mx_cv_s["qaly"] <- qaly * alive_midyear
    mx_cv_s["qaly"] <- qaly * (1 + d_midyear)
    
    ### update mx_cv ###############################
    
    mx_cv[s, ] <- mx_cv_s
    
  }
  
  ### update v with cumulative probabilities  ####################
  
  vP_1_temp <- vP_1[mx_cv_rownames]
  # only update outcomes that are dealt with in mx_cv
  for (outcome in mx_cv_colnames){
    vP_2[outcome] <-  vP_1_temp %*% mx_cv[, outcome]
  }
  
  vP_2["alive_at_end"] <- vP_2["alive_at_start"] - vP_2["endpt_VD"] - vP_2["endpt_NVD"]
  
  ### return ######################################
  
  retval <- list(vP_2 = vP_2)
  return(retval) 
}

###################################################################
### wrappers
###################################################################

master <- function( 
  data_for_extrapolation_dir, data_risk_equations_dir,
  output_dir,
  states_info, endpts_CV, names_mx_mon,
  pop_suffix, secondary,
  tp_survival_nvd, tp_survival_vd,
  tp_survival_vd_stroke, tp_survival_vd_stroke_mi, tp_survival_hf,
  coeff_nvd, anc_nvd,
  coeff_vd, anc_vd,
  coeff_vd_stroke, anc_vd_stroke,
  coeff_vd_stroke_mi, anc_vd_stroke_mi,
  coeff_hf, anc_hf,
  tx_statin, tx_antihyp, tx_antip,
  mx_ckd, alb, Q, Q_ckd, tp_esrd, covmeans,
  coeff_qol,
  mx_b, mx_t, p_b,
  df_tx_flag, ckd_clin, hyp_bp, update_diab,
  stop_expr,
  statin_diab, statin_nvd, statin_qol_decr, aki_incidence, aki_qol_decr,
  j_G2, j_G3a, j_G3b, j_G4, j_G5, 
  n_cores) {
  
  # beta * X for baseline characteristics
  vars_b <- setdiff(colnames(mx_b), c("patid_temp", "ckd_base"))
  lamRE <- prelimRiskEquations(vars_b = vars_b,
                               mx_b = mx_b,
                               coeff_NVD = coeff_nvd,
                               coeff_VD = coeff_vd, coeff_VD_Stroke = coeff_vd_stroke,
                               coeff_MI_Stroke_VD = coeff_vd_stroke_mi,
                               coeff_HF = coeff_hf)
  
  
  lamRE_lambdas <- lamRE$lambdas
  lamRE_coeff <- lamRE$coeff
  
  coeff_NVD_t <- lamRE_coeff$coeff_NVD$coeff_t
  coeff_VD_t <- lamRE_coeff$coeff_VD$coeff_t
  coeff_VD_Stroke_t <- lamRE_coeff$coeff_VD_Stroke$coeff_t
  coeff_MI_Stroke_VD_t <- lamRE_coeff$coeff_MI_Stroke_VD$coeff_t
  coeff_HF_t <- lamRE_coeff$coeff_HF$coeff_t
  
  # Preliminary Qol
  qolRE <- prelimQol(vars_b = vars_b, mx_b = mx_b, coeff_qol = coeff_qol)
  qolRE_lambdas <- qolRE$lambdas
  qolRE_coeff <- qolRE$coeff
  coeff_qol_t <- qolRE_coeff$coeff_qol$coeff_t
  
  # read off names
  coeff_names_cost_qol <- setdiff(
    #union(names(coeff_qol_t), names(coeff_cost_t)),
    names(coeff_qol_t),
    c("age_cent", "diab1", "diab2", "ckd3a", "ckd3b", "ckd4", "ckd5", "dialysis", "transplant"))
  coeff_names_cost_qol <- setdiff(coeff_names_cost_qol, "prev_HF")
  
  ### identify patients
  
  ids <- 1 : nrow(mx_b)
  #ids <- 1 : 10
  # set.seed(240485)
  # ids <- sort(sample(1 : nrow(mx_b), 10))
  # 
  #i <- ids[1]
  
  ### start a cluster
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("H_g", "H_w", "gen_nextCycle",
                      "get_pcond", "tp_survival_g", "tp_survival_w",
                      ".lamVec", ".lamM", "get_tp_mx_ckd", #"cycle_cost",
                      "cycle_qol", "get_state_combined")) 
  clusterEvalQ(cl, library(expm))
  registerDoParallel(cl)
  
  #retval <- list()
  #for (i in ids){
  #print(i)
  retval <- foreach (i = ids) %dopar% {
    
    ### extract information on the patients
    
    # baseline characteristics
    vX_b <- mx_b[i, ]
    vX_i <- mx_t[i, ]
    vP_b <- p_b[i, ]
    treat_i <- df_tx_flag[df_tx_flag[, "patid_temp"] == vX_b["patid_temp"], ]
    ckd_clin_i <- ckd_clin[i, "ckd_clin"]
    hyp_bp_i <- hyp_bp[i, "hyp_bp"]
    update_diab_i <- update_diab[i, "update_diab"]
    
    # baseline beta * X for risk equations
    lam_NVD_b <- lamRE_lambdas$lam_NVD_b_All[i]
    lam_VD_b <- lamRE_lambdas$lam_VD_b_All[i]
    lam_VD_Stroke_b <- lamRE_lambdas$lam_VD_Stroke_b_All[i]
    lam_MI_Stroke_VD_b <- lamRE_lambdas$lam_MI_Stroke_VD_b_All[i]
    lam_HF_b <- lamRE_lambdas$lam_HF_b_All[i]
    lam_qol_b <- qolRE_lambdas$lam_qol_b[i]
    
    # ckd model
    alb_i <- alb[i, "acr"] + 1
    covmeans_i <- covmeans[[alb_i]]
    Q_ckd_i <- Q_ckd[[alb_i]]
    mx_ckd_i <- mx_ckd[i, ]
    
    # generate output matrix
    N <- max(eval(stop_expr) - 1, 1)
    output <- matrix(nrow = N, ncol = length(vP_b))
    
    #j <- 1
    
    # perform the simulation
    for (j in 1 : N) {
      
      #print(j)
      
      alpha <- gen_nextCycle(j = j, secondary = secondary, states_info = states_info, endpts_CV = endpts_CV, names_mx_mon = names_mx_mon,
                             tp_survival_nvd = tp_survival_nvd, coeff_NVD_t = coeff_NVD_t, anc_nvd = anc_nvd, lam_NVD_b = lam_NVD_b,
                             tp_survival_vd = tp_survival_vd, coeff_VD_t = coeff_VD_t, anc_vd = anc_vd, lam_VD_b = lam_VD_b,
                             tp_survival_vd_stroke = tp_survival_vd_stroke, coeff_VD_Stroke_t = coeff_VD_Stroke_t, anc_vd_stroke = anc_vd_stroke, lam_VD_Stroke_b = lam_VD_Stroke_b,
                             tp_survival_vd_stroke_mi = tp_survival_vd_stroke_mi, coeff_MI_Stroke_VD_t = coeff_MI_Stroke_VD_t,
                             anc_vd_stroke_mi = anc_vd_stroke_mi, lam_MI_Stroke_VD_b = lam_MI_Stroke_VD_b,
                             tp_survival_hf = tp_survival_hf, coeff_HF_t = coeff_HF_t, anc_hf = anc_hf, lam_HF_b = lam_HF_b,
                             tx_statin = tx_statin, tx_antihyp = tx_antihyp, tx_antip = tx_antip,
                             coeff_qol_t = coeff_qol_t, lam_qol_b = lam_qol_b,
                             coeff_names_cost_qol = coeff_names_cost_qol,
                             vP_b = vP_b, vX_i = vX_i, 
                             covmeans_i = covmeans_i, Q_ckd_i = Q_ckd_i, mx_ckd_i = mx_ckd_i, tp_esrd = tp_esrd, Q = Q,
                             treat_i = treat_i, ckd_clin_i = ckd_clin_i, hyp_bp_i = hyp_bp_i, update_diab_i = update_diab_i,
                             j_G2 = j_G2, j_G3a = j_G3a, j_G3b = j_G3b, j_G4 = j_G4, j_G5 = j_G5, 
                             statin_diab = statin_diab, statin_nvd = statin_nvd, statin_qol_decr = statin_qol_decr,
                             aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr)
      
      vP_b <- alpha[["vP_2"]]
      output[j, ] <- vP_b
    }
    # save the output
    #retval[[i]] <- output
    return(output)
  }
  stopCluster(cl)
  
  ### combine output
  
  df_output <- as.data.frame(do.call("rbind", retval))
  colnames(df_output) <- colnames(p_b)
  
  ####### cleanup the dataset #######
  
  # non-fatal combined endpoints
  df_output <- df_output %>%
    arrange(patid_temp, cycle) %>%
    group_by(patid_temp) %>%
    mutate(ly = alive_at_start - 1/2 * (endpt_NVD + endpt_VD)) %>%
    select(patid_temp, ly, qaly) %>% 
    summarise(ly = sum(ly), qaly = sum(qaly))
  
  return(df_output) 
}

wrapper_est <- function(
  data_for_extrapolation_dir,
  data_risk_equations_dir,
  output_dir,
  mx_b_prefix, mx_t_prefix, p_b_prefix, state_prefix, 
  mx_ckd_prefix, alb_prefix, Q_prefix,
  df_tx_flag_prefix, ckd_clin_prefix, hyp_bp_prefix, update_diab_prefix,
  stop_expr,
  statin_diab, statin_nvd, statin_qol_decr, aki_incidence, aki_qol_decr,
  output_file_prefix,
  cf_prefix, tx, n_cores
) {
  
  # states & events
  setwd(data_for_extrapolation_dir)
  states <- get(load(paste0(states_prefix, ".Rdata")))
  states_info <- states$states_info
  endpts_CV <- states$endpts_CV
  # remove irrelevant endpoints
  # perhaps uglier than saving a separate file but only needs to be done once
  endpts_CV <- setdiff(endpts_CV, grep("not", endpts_CV, value = TRUE))
  
  # other parameters
  covmeans <- get(load(paste0(covmeans_prefix, ".Rdata")))
  Q <- get(load(paste0(Q_prefix, ".Rdata")))
  
  # survival distributions
  tp_survival_nvd <- "tp_survival_g"
  tp_survival_vd <- "tp_survival_g" 
  tp_survival_vd_stroke <- "tp_survival_w"
  tp_survival_vd_stroke_mi <- "tp_survival_w"
  tp_survival_hf <- "tp_survival_w"
  
  # read off indices for the CKD stages
  # ugly as need to load the file that's not needed but the alternative is to load it 4,000 times
  p_temp <- get(load(paste0(p_b_prefix, "F_0.Rdata")))
  j_1_1 <- which(colnames(p_temp) == "1_1")
  j_G2 <- j_1_1:(j_1_1 + 6)
  j_G3a <- j_G2 + 7
  j_G3b <- j_G3a + 7
  j_G4 <- j_G3b + 7
  j_G5 <- j_G4 + 7
  rm(p_temp)
  
  # other parameters that only need to be defined once
  names_mx_mon <- as.vector(outer(c(1:7), c(1:5), function(x, y) paste(y, x, sep = "_")))
  
  for (sex in c("F", "M")) {
    
    for (secondary in c(0, 1)){
      
      pop_suffix <- paste(sex, secondary, sep = "_")
      # coefficients
      setwd(file.path(data_risk_equations_dir))
      cf_temp <- get(load(paste0(cf_prefix, "_", pop_suffix, ".Rdata")))
      # load patients characteristics & states
      setwd(data_for_extrapolation_dir)
      mx_b <- get(load(paste0(mx_b_prefix, pop_suffix, ".Rdata")))
      mx_t <- get(load(paste0(mx_t_prefix, pop_suffix, ".Rdata")))
      p_b <- get(load(paste0(p_b_prefix, pop_suffix, ".Rdata")))
      mx_ckd <- get(load(paste0(mx_ckd_prefix, pop_suffix, ".Rdata")))
      alb <- get(load(paste0(alb_prefix, pop_suffix, ".Rdata")))
      df_tx_flag <- get(load(paste0(df_tx_flag_prefix, pop_suffix, ".Rdata")))
      ckd_clin <- get(load(paste0(ckd_clin_prefix, pop_suffix, ".Rdata")))
      hyp_bp <- get(load(paste0(hyp_bp_prefix, pop_suffix, ".Rdata")))
      update_diab <- get(load(paste0(update_diab_prefix, pop_suffix, ".Rdata")))
      
      # nonvascular death
      temp <- cf_temp[["nvd"]]
      coeff_nvd <- temp[[1]]
      anc_nvd <- temp[[2]]
      
      # vascular death
      temp <- cf_temp[["vd"]]
      coeff_vd <- temp[[1]]
      anc_vd <- temp[[2]]
      
      # stroke or vd
      temp <- cf_temp[["vd_stroke"]]
      coeff_vd_stroke <- temp[[1]]
      anc_vd_stroke <- temp[[2]]
      
      # mi, stroke or vd
      temp <- cf_temp[["vd_stroke_mi"]]
      coeff_vd_stroke_mi <- temp[[1]]
      anc_vd_stroke_mi <- temp[[2]]
      
      # heart failure
      temp <- cf_temp[["hf"]]
      coeff_hf <- temp[[1]]
      anc_hf <- temp[[2]]
      
      # ckd model
      Q_ckd <- list(cf_temp[["ckd_unmeasured"]], cf_temp[["ckd_normal"]], cf_temp[["ckd_micro"]], cf_temp[["ckd_macro"]])
      tp_esrd <- cf_temp[["tp_esrd"]]
      
      # qol
      coeff_qol <- cf_temp[["qol"]]
      
      # treatment effects
      
      tx_statin_vd_a20 <- tx[["tx_statin_vd_a20"]]
      tx_statin_stroke_a20 <- tx[["tx_statin_stroke_a20"]]
      tx_statin_mi_a20 <- tx[["tx_statin_mi_a20"]]
      tx_statin_1 <- rbind(tx_statin_vd_a20, tx_statin_stroke_a20, tx_statin_mi_a20)
      
      tx_statin_vd_a80 <- tx[["tx_statin_vd_a80"]]
      tx_statin_stroke_a80 <- tx[["tx_statin_stroke_a80"]]
      tx_statin_mi_a80 <- tx[["tx_statin_mi_a80"]]
      tx_statin_2 <- rbind(tx_statin_vd_a80, tx_statin_stroke_a80, tx_statin_mi_a80)
      
      tx_statin <- list(prim = tx_statin_1, sec = tx_statin_2)
      
      tx_antihyp_vd <- tx[["tx_antihyp_vd"]]
      tx_antihyp_stroke <- tx[["tx_antihyp_stroke"]]
      tx_antihyp_mi <- tx[["tx_antihyp_mi"]]
      tx_antihyp_hf <- tx[["tx_antihyp_hf"]]
      tx_antihyp <- rbind(tx_antihyp_vd, tx_antihyp_stroke, tx_antihyp_mi, tx_antihyp_hf)
      
      tx_antip_vd <- tx[["tx_antip_vd"]]
      tx_antip_stroke <- tx[["tx_antip_stroke"]]
      tx_antip_mi <- tx[["tx_antip_mi"]]
      tx_antip <- rbind(tx_antip_vd, tx_antip_stroke, tx_antip_mi)
      
      ### run the master function
      alpha <- master( 
        data_for_extrapolation_dir = data_for_extrapolation_dir,
        data_risk_equations_dir = data_risk_equations_dir,
        output_dir = output_dir,
        states_info = states_info, endpts_CV = endpts_CV,
        stop_expr = stop_expr, names_mx_mon = names_mx_mon,
        tp_survival_nvd = tp_survival_nvd, tp_survival_vd = tp_survival_vd,
        tp_survival_vd_stroke = tp_survival_vd_stroke, tp_survival_vd_stroke_mi = tp_survival_vd_stroke_mi,
        tp_survival_hf = tp_survival_hf,
        coeff_nvd = coeff_nvd, anc_nvd = anc_nvd,
        coeff_vd = coeff_vd, anc_vd = anc_vd,
        coeff_vd_stroke = coeff_vd_stroke, anc_vd_stroke = anc_vd_stroke,
        coeff_vd_stroke_mi = coeff_vd_stroke_mi, anc_vd_stroke_mi = anc_vd_stroke_mi,
        coeff_hf = coeff_hf, anc_hf = anc_hf,
        tx_statin = tx_statin, tx_antihyp = tx_antihyp, tx_antip = tx_antip,
        Q_ckd = Q_ckd, tp_esrd = tp_esrd, covmeans = covmeans,
        coeff_qol = coeff_qol,
        pop_suffix = pop_suffix, secondary = secondary,
        mx_b = mx_b, mx_t = mx_t, p_b = p_b, mx_ckd = mx_ckd, alb = alb, Q = Q,
        df_tx_flag = df_tx_flag, ckd_clin = ckd_clin, hyp_bp = hyp_bp, update_diab = update_diab,
        statin_diab = statin_diab, statin_nvd = statin_nvd, statin_qol_decr = statin_qol_decr,
        j_G2 = j_G2, j_G3a = j_G3a, j_G3b = j_G3b, j_G4 = j_G4, j_G5 = j_G5, 
        aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr,
        n_cores = n_cores)
      # save the output
      setwd(output_dir)
      save(alpha, file = paste0(output_file_prefix, "_", pop_suffix, ".Rdata"))
    }
  }
}

wrapper_psa <- function(
  data_for_extrapolation_dir,
  data_risk_equations_dir,
  output_dir,
  mx_b_prefix, mx_t_prefix, p_b_prefix, state_prefix, 
  mx_ckd_prefix, alb_prefix, Q_prefix,
  df_tx_flag_prefix, ckd_clin_prefix, hyp_bp_prefix, update_diab_prefix,
  stop_expr,
  statin_diab, statin_nvd, statin_qol_decr, aki_incidence, aki_qol_decr,
  output_file_prefix,
  boot_prefix, n_sim_start, n_sim_end, n_cores
) {
  
  # states & events
  setwd(data_for_extrapolation_dir)
  states <- get(load(paste0(states_prefix, ".Rdata")))
  states_info <- states$states_info
  endpts_CV <- states$endpts_CV
  # remove irrelevant endpoints
  # perhaps uglier than saving a separate file but only needs to be done once
  endpts_CV <- setdiff(endpts_CV, grep("not", endpts_CV, value = TRUE))
  
  # other parameters
  covmeans <- get(load(paste0(covmeans_prefix, ".Rdata")))
  Q <- get(load(paste0(Q_prefix, ".Rdata")))
  
  # survival distributions
  tp_survival_nvd <- "tp_survival_g"
  tp_survival_vd <- "tp_survival_g" 
  tp_survival_vd_stroke <- "tp_survival_w"
  tp_survival_vd_stroke_mi <- "tp_survival_w"
  tp_survival_hf <- "tp_survival_w"
  
  # read off indices for the CKD stages
  # ugly as need to load the file that's not needed but the alternative is to load it 4,000 times
  p_temp <- get(load(paste0(p_b_prefix, "F_0.Rdata")))
  j_1_1 <- which(colnames(p_temp) == "1_1")
  j_G2 <- j_1_1:(j_1_1 + 6)
  j_G3a <- j_G2 + 7
  j_G3b <- j_G3a + 7
  j_G4 <- j_G3b + 7
  j_G5 <- j_G4 + 7
  rm(p_temp)
  
  # other parameters that only need to be defined once
  names_mx_mon <- as.vector(outer(c(1:7), c(1:5), function(x, y) paste(y, x, sep = "_")))
  
  for (sex in c("F", "M")) {
    
    for (secondary in c(0, 1)){
      
      pop_suffix <- paste(sex, secondary, sep = "_")
      # coefficients
      setwd(file.path(data_risk_equations_dir, "psa_input"))
      cf_boot_temp <- get(load(paste0(boot_prefix, "_", pop_suffix, ".Rdata")))
      # load patients characteristics & states
      setwd(data_for_extrapolation_dir)
      mx_b <- get(load(paste0(mx_b_prefix, pop_suffix, ".Rdata")))
      mx_t <- get(load(paste0(mx_t_prefix, pop_suffix, ".Rdata")))
      p_b <- get(load(paste0(p_b_prefix, pop_suffix, ".Rdata")))
      mx_ckd <- get(load(paste0(mx_ckd_prefix, pop_suffix, ".Rdata")))
      alb <- get(load(paste0(alb_prefix, pop_suffix, ".Rdata")))
      df_tx_flag <- get(load(paste0(df_tx_flag_prefix, pop_suffix, ".Rdata")))
      ckd_clin <- get(load(paste0(ckd_clin_prefix, pop_suffix, ".Rdata")))
      hyp_bp <- get(load(paste0(hyp_bp_prefix, pop_suffix, ".Rdata")))
      update_diab <- get(load(paste0(update_diab_prefix, pop_suffix, ".Rdata")))
      
      for (n in n_sim_start : n_sim_end){
        
        print(n)
        
        cf_boot_n <- cf_boot_temp[[n]]
        
        # nonvascualr death
        temp <- cf_boot_n[["nvd"]]
        coeff_nvd <- temp[[1]]
        anc_nvd <- temp[[2]]
        
        # vascular death
        temp <- cf_boot_n[["vd"]]
        coeff_vd <- temp[[1]]
        anc_vd <- temp[[2]]
        
        # stroke or vd
        temp <- cf_boot_n[["vd_stroke"]]
        coeff_vd_stroke <- temp[[1]]
        anc_vd_stroke <- temp[[2]]
        
        # mi, stroke or vd
        temp <- cf_boot_n[["vd_stroke_mi"]]
        coeff_vd_stroke_mi <- temp[[1]]
        anc_vd_stroke_mi <- temp[[2]]
        
        # heart failure
        temp <- cf_boot_n[["hf"]]
        coeff_hf <- temp[[1]]
        anc_hf <- temp[[2]]
        
        # ckd model
        Q_ckd <- list(cf_boot_n[["ckd_unmeasured"]], cf_boot_n[["ckd_normal"]], cf_boot_n[["ckd_micro"]], cf_boot_n[["ckd_macro"]])
        tp_esrd <- cf_boot_n[["tp_esrd"]]
        
        # qol
        coeff_qol <- cf_boot_n[["qol"]]
        
        # treatment effects
        
        tx_statin_vd_a20 <- cf_boot_n[["tx_statin_vd_a20"]]
        tx_statin_stroke_a20 <- cf_boot_n[["tx_statin_stroke_a20"]]
        tx_statin_mi_a20 <- cf_boot_n[["tx_statin_mi_a20"]]
        tx_statin_1 <- rbind(tx_statin_vd_a20, tx_statin_stroke_a20, tx_statin_mi_a20)
        
        tx_statin_vd_a80 <- cf_boot_n[["tx_statin_vd_a80"]]
        tx_statin_stroke_a80 <- cf_boot_n[["tx_statin_stroke_a80"]]
        tx_statin_mi_a80 <- cf_boot_n[["tx_statin_mi_a80"]]
        tx_statin_2 <- rbind(tx_statin_vd_a80, tx_statin_stroke_a80, tx_statin_mi_a80)
        
        tx_statin <- list(prim = tx_statin_1, sec = tx_statin_2)
        
        tx_antihyp_vd <- cf_boot_n[["tx_antihyp_vd"]]
        tx_antihyp_stroke <- cf_boot_n[["tx_antihyp_stroke"]]
        tx_antihyp_mi <- cf_boot_n[["tx_antihyp_mi"]]
        tx_antihyp_hf <- cf_boot_n[["tx_antihyp_hf"]]
        tx_antihyp <- rbind(tx_antihyp_vd, tx_antihyp_stroke, tx_antihyp_mi, tx_antihyp_hf)
        
        tx_antip_vd <- cf_boot_n[["tx_antip_vd"]]
        tx_antip_stroke <- cf_boot_n[["tx_antip_stroke"]]
        tx_antip_mi <- cf_boot_n[["tx_antip_mi"]]
        tx_antip <- rbind(tx_antip_vd, tx_antip_stroke, tx_antip_mi)
        
        ### run the master function
        alpha <- master( 
          data_for_extrapolation_dir = data_for_extrapolation_dir,
          data_risk_equations_dir = data_risk_equations_dir,
          output_dir = output_dir,
          states_info = states_info, endpts_CV = endpts_CV,
          stop_expr = stop_expr, names_mx_mon = names_mx_mon,
          tp_survival_nvd = tp_survival_nvd, tp_survival_vd = tp_survival_vd,
          tp_survival_vd_stroke = tp_survival_vd_stroke, tp_survival_vd_stroke_mi = tp_survival_vd_stroke_mi,
          tp_survival_hf = tp_survival_hf,
          coeff_nvd = coeff_nvd, anc_nvd = anc_nvd,
          coeff_vd = coeff_vd, anc_vd = anc_vd,
          coeff_vd_stroke = coeff_vd_stroke, anc_vd_stroke = anc_vd_stroke,
          coeff_vd_stroke_mi = coeff_vd_stroke_mi, anc_vd_stroke_mi = anc_vd_stroke_mi,
          coeff_hf = coeff_hf, anc_hf = anc_hf,
          tx_statin = tx_statin, tx_antihyp = tx_antihyp, tx_antip = tx_antip,
          Q_ckd = Q_ckd, tp_esrd = tp_esrd, covmeans = covmeans,
          coeff_qol = coeff_qol,
          pop_suffix = pop_suffix, secondary = secondary,
          mx_b = mx_b, mx_t = mx_t, p_b = p_b, mx_ckd = mx_ckd, alb = alb, Q = Q,
          df_tx_flag = df_tx_flag, ckd_clin = ckd_clin, hyp_bp = hyp_bp, update_diab = update_diab,
          statin_diab = statin_diab, statin_nvd = statin_nvd, statin_qol_decr = statin_qol_decr,
          j_G2 = j_G2, j_G3a = j_G3a, j_G3b = j_G3b, j_G4 = j_G4, j_G5 = j_G5, 
          aki_incidence = aki_incidence, aki_qol_decr = aki_qol_decr,
          n_cores = n_cores)
        # save the output
        setwd(output_dir)
        save(alpha, file = paste0(output_file_prefix, "_", pop_suffix, "_", n, ".Rdata"))
      }
    }
  }
}
