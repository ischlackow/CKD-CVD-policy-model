###################################################################
#   Harrell's c-statistic
###################################################################

###########################################
### Prepare the data
###########################################

rm(list = ls())

library(gdata)
library(reshape2)
library(ggplot2)
library(scales)
library(xtable)
library(plyr)
library(survival)
library(data.table)

source(file.path("K:", "SHARP_EE", "Monitoring CKD CHF", "ckd-monitoring", "src", "model_paper", "model_paper_functions.R"))

sourceDirPred <- "K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\output\\lifetime_datasets\\"

### actual data
df <- data.table(get(load("K:/SHARP_EE/Monitoring CKD CHF/Data/Derived CPRD data/df_Aevents_cross_validation.Rdata")))
dfAct <- df

# dataset
alpha <- subset(dfAct, select = c("patid", "estimation", "ckd", 
                                  "NVD", "TTE_NVD",
                                  "VD", "TTE_VD",
                                  "TTE_HF_hosp", "HF_hosp",
                                  "TTE_Stroke_or_VD", "Stroke_or_VD",
                                  "TTE_MI_Stroke_or_VD", "MI_Stroke_or_VD"))

alpha <- alpha[,lapply(.SD, sum), .SDcols=c("NVD", "VD", "Stroke_or_VD", "MI_Stroke_or_VD", "HF_hosp", "TTE_NVD", "TTE_VD", "TTE_MI_Stroke_or_VD", 
                                                "TTE_Stroke_or_VD", "TTE_HF_hosp"), by=c("patid", "ckd", "estimation")]
#alpha <- alpha[estimation==1,]
alpha <- alpha[estimation==0,]
alpha <- data.frame(alpha)
# create survival objects
#alpha$VD[alpha$FU_VD > 5] <- 0
#alpha$FU_VD[alpha$FU_VD > 5] <- 5
alpha$surv_VD <- (Surv(alpha$TTE_VD, alpha$VD))
alpha$surv_NVD <- (Surv(alpha$TTE_NVD, alpha$NVD))

alpha$surv_Stroke_or_VD <- (Surv(alpha$TTE_Stroke_or_VD, alpha$Stroke_or_VD))
alpha$surv_MI_Stroke_or_VD <- (Surv(alpha$TTE_MI_Stroke_or_VD, alpha$MI_Stroke_or_VD))

alpha$surv_HF_hosp <- (Surv(alpha$TTE_HF_hosp, alpha$HF_hosp))

head(alpha)

alpha1 <- data.frame(alpha) #data.table(alpha)
ids <- alpha1$patid
rm(alpha)
rm(df, df_Aevents, dfAct)

### predicted data

setwd(sourceDirPred)
dfPred <- data.table(.combine_4_into_1(file_prefix = "df_output_int_val", file_suffix = NULL, dir = sourceDirPred))

df_id_all <- data.table(get(load("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\data\\processed\\df_id_all.Rdata")))
head(df_id_all)
df_id_all <- df_id_all[,.(patid, patid_temp)]
dfPred <- merge(dfPred, df_id_all, by="patid_temp", all.x=TRUE)
dfPred <- dfPred[patid %in% ids,]

### Temp for baseline data currently missing:
load("K:/SHARP_EE/Monitoring CKD CHF/ckd-monitoring/data/for_extrapolation/mx_b_int_val_F_0.Rdata")
load("K:/SHARP_EE/Monitoring CKD CHF/ckd-monitoring/data/for_extrapolation/mx_b_int_val_F_1.Rdata")
load("K:/SHARP_EE/Monitoring CKD CHF/ckd-monitoring/data/for_extrapolation/mx_b_int_val_M_0.Rdata")
load("K:/SHARP_EE/Monitoring CKD CHF/ckd-monitoring/data/for_extrapolation/mx_b_int_val_M_1.Rdata")

baseline <- data.table(rbind(mx_b_int_val_F_0, mx_b_int_val_F_1, mx_b_int_val_M_0, mx_b_int_val_M_1))
names(baseline)
baseline <- baseline[,.(patid_temp, ckd_base)]

dfPred <- merge(dfPred, baseline, by="patid_temp", all.x = TRUE)

rm(baseline, mx_b_int_val_F_0, mx_b_int_val_F_1, mx_b_int_val_M_0, mx_b_int_val_M_1)

dfPred <- subset(dfPred, select = c("patid", "cycle","ckd_base", "prev_HF", "HF_hosp",
                                    grep("endpt_", colnames(dfPred), value = TRUE, ignore.case = TRUE),
                                    grep("state_", colnames(dfPred), value = TRUE, ignore.case = TRUE),
                                    grep("alive", colnames(dfPred), value = TRUE, ignore.case = TRUE),
                                    grep("1", colnames(dfPred), value = TRUE, ignore.case = TRUE),
                                    grep("2", colnames(dfPred), value = TRUE, ignore.case = TRUE),
                                    grep("4", colnames(dfPred), value = TRUE, ignore.case = TRUE),
                                    grep("fu", colnames(dfPred), value = TRUE, ignore.case = TRUE)))
head(dfPred)

addCols <- function(colNames) {
  # check for another way?
  return(eval(parse(text = paste("dfPred[, '", colNames, "']", sep = "", collapse = "+"))))
}

dfPred[,FU_NVD:= dfPred$alive_at_start]
dfPred[,FU_VD := alive_at_start - endpt_NVD]
dfPred[,FU_S_VD_temp := fu_s_vd_temp]
dfPred[,FU_MI_S_VD_temp := fu_mi_s_vd_temp]
dfPred[,FU_HF := alive_at_start - prev_HF]

dfPred[,FU_S_VD:=shift(FU_S_VD_temp,1), by=.(patid)]
dfPred[is.na(FU_S_VD), FU_S_VD:=1]

dfPred[,FU_MI_S_VD:=shift(FU_MI_S_VD_temp,1), by=.(patid)]
dfPred[is.na(FU_MI_S_VD), FU_MI_S_VD:=1]


dfPred[,ckd := ifelse(dfPred$ckd_base==1, "ckd2",
                      ifelse(dfPred$ckd_base==2, "ckd3a", 
                             ifelse(dfPred$ckd_base==3, "ckd3b", "ckd4/5"))) ]
dfPred[,ckd := factor(dfPred$ckd, labels=c("ckd2", "ckd3a", "ckd3b", "ckd4/5"))]

# calculate KM-type event rates

sort(colnames(dfPred))

beta <- subset(dfPred, cycle >= 1, select = c("patid", "cycle", "ckd",
                                              "endpt_NVD", "FU_NVD",
                                              "endpt_VD", "FU_VD",
                                              "endpt_vd_stroke", "FU_S_VD",
                                              "endpt_vd_stroke_mi", "FU_MI_S_VD",
                                              "HF_hosp", "FU_HF"))

# rates
beta1 <- beta
beta1$NVD <- 1 - beta1$endpt_NVD / beta1$FU_NVD
beta1$VD <- 1 - beta1$endpt_VD / beta1$FU_VD
beta1$Stroke_VD <- 1 - beta1$endpt_vd_stroke/beta1$FU_S_VD
beta1$MI_Stroke_VD <- 1 - beta1$endpt_vd_stroke_mi/beta1$FU_MI_S_VD
beta1$HF_hosp <- 1 - beta1$HF_hosp/beta1$FU_HF

rm(beta, dfPred)

## KM product
beta1 <- data.table(beta1)
beta1 <-  beta1[patid %in% ids]

beta1[,`:=`(VD_pred = cumprod(VD),    
             NVD_pred = cumprod(NVD),
             HF_pred = cumprod(HF_hosp),
             Stroke_VD_pred = cumprod(Stroke_VD),
             MI_Stroke_VD_pred = cumprod(MI_Stroke_VD)), by=patid]

# clean
beta3 <- beta1[cycle==5, .(patid, NVD_pred, VD_pred, HF_pred, Stroke_VD_pred, MI_Stroke_VD_pred)]
beta3 <- data.frame(beta3)

rm(beta1)
# combine
alpha1 <- data.frame(alpha1)
gamma <- merge(alpha1, beta3, by = "patid", all.y=TRUE)
gamma

############################
### calculate c-statistic
############################

library(Hmisc)
library(xlsx)# 

endpts <- list(
  list(lab = "Vascular death", colname_act = "surv_VD", colname_pred = "VD_pred"),
  list(lab = "Non-Vascular death", colname_act = "surv_NVD", colname_pred = "NVD_pred"),
  list(lab = "Stroke or vascular death", colname_act = "surv_Stroke_or_VD", colname_pred = "Stroke_VD_pred"),
  list(lab = "MI, Stroke or Vascular Death", colname_act = "surv_MI_Stroke_or_VD", colname_pred = "MI_Stroke_VD_pred"),
  list(lab = "HF hospitalisation", colname_act = "surv_HF_hosp", colname_pred = "HF_pred")
)

# gamma <- output
output <-   data.frame(subgroup = character(), endpt = character(), cstat = numeric(), cstat_L = numeric(), cstat_U = numeric())
rm(alpha1, beta, beta1, beta3, dfPred, gamma1)
# by CKD stage
subgroups <- c("ckd2", "ckd3a", "ckd3b", "ckd4", "ckd5")
subgroup <- subgroups[1]
endpt <- endpts[[1]]

for (subgroup in subgroups){
  print(subgroup)
  for (endpt in endpts) {
    print(endpt)
    temp <- if(subgroup=="all") gamma else
      subset(gamma, ckd == subgroup)
    x <- rcorr.cens(temp[, endpt$colname_pred], temp[, endpt$colname_act])
    cstat <- x["C Index"]
    se_const <- 1.96 * x["S.D."] / 2
    output <- rbind(output, data.frame(subgroup = subgroup, endpt = endpt$lab, 
                                       cstat = cstat, cstat_L = cstat - se_const, cstat_U = cstat + se_const))
  }
}
output
test <- melt(output, id.vars = c("subgroup", "endpt"))
test <- subset(test, endpt %in% c("Vascular death", "Stroke or vascular death", "MI, Stroke or Vascular Death"))
test2 <- dcast(test, subgroup ~ endpt + variable)
#write.csv(test2, file="K:\\SHARP_EE\\Monitoring CKD CHF\\Documents\\Papers\\Model Paper\\outputs\\C_Statistic_estimation.csv", row.names = FALSE)
write.csv(test2, file="K:\\SHARP_EE\\Monitoring CKD CHF\\Documents\\Papers\\Model Paper\\outputs\\C_Statistic_validation.csv", row.names = FALSE)