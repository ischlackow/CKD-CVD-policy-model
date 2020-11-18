#############################################################
### Random numbers for the manuscript
#############################################################

rm(list = ls())

library(tidyverse)
output_dir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Documents\\Papers\\Model Paper\\outputs\\"


df <- get(load("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\data\\processed\\df_Aevents_cross_validation.Rdata"))

# baseline info
df <- df %>% 
  filter(cycle == 0) %>%
  mutate(baseline_CV = pmin(Base_Cardiovascular + Base_Cerebrovascular + Base_HF, 1),
         ckd_clin = as.numeric(ckd != "ckd2" | acr > 1 | crd == 1))

# total number
nrow(df)
round(nrow(df) / 1000000, digits = 1)

# male / female distribution
df %>%
  group_by(sex) %>%
  summarise(n = round(100 * n() / nrow(df), digits = 0))

# secondary / primary distribution
df %>%
  group_by(baseline_CV) %>%
  summarise(n = round(100 * n() / nrow(df), digits = 0))

# eGFR stage distribution
df %>%
  group_by(ckd_updated) %>%
  summarise(n = round(100 * n() / nrow(df), digits = 0))
df %>%
  group_by(ckd_updated) %>%
  summarise(n = round(100 * n() / nrow(df), digits = 1))
df %>%
  group_by(ckd) %>%
  summarise(n = round(100 * n() / nrow(df), digits = 0))

# follow-up
mean(df$FU)
round(summary(df$FU), digits = 1)

# number of practices
length(unique(df$pracid))

# missing data
sort(colnames(df))
for (v in  c("ethnicity", "smoking", "bmi", "imd", "sbp", "dbp", "chol_hdl", "acr")) {
  print(v)
  print(sum(is.na(df[, v])))
  print(round(100 * sum(is.na(df[, v])) / nrow(df), digits = 1))
}

# missing data by eGFR stage
sort(colnames(df))


v <- c("ethnicity", "smoking", "bmi", "imd", "sbp", "dbp", "chol_hdl", "acr_2")
df_na <- df %>%
  mutate(acr_2 = case_when(acr == 0 ~ NA,
                           TRUE ~ TRUE)) %>%
  select_at(c("ckd", v)) %>%
  group_by(ckd) %>%
  summarise_all( ~ str_c(scales::number(sum(is.na(.x)), big.mark = ","), " (", round(100 * sum(is.na(.x)) / n(), 2), "%)"))

df_na

df_na_t <- df_na %>% rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  pivot_wider(names_from = rowname, values_from = value)
setwd(output_dir)
write_csv(df_na_t, path = "missing_by_eGFR.csv")

# acr variable
df %>%
  group_by(acr) %>%
  summarise(n = n(), p = round(100 * n() / nrow(df), digits = 1))

# age distribution
df %>%
  group_by(ckd) %>%
  summarise(age = mean(age))

#############################################################
### Demographics
#############################################################

rm(list = ls())

library(tidyverse)
library(data.table)

source("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\src\\model_paper\\numbers_functions.R")

outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Documents\\Papers\\Model Paper\\outputs\\"

df <- data.table(get(load("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\data\\processed\\df_Aevents_cross_validation.Rdata")))
df <- df[cycle == 0]

df2 <- data.table(get(load("K:/SHARP_EE/Monitoring CKD CHF/Data/Derived CPRD data/df_incBase_hist_incACR.Rdata")))
table(df2$ckd_READ, exclude = NULL)
nrow(df)
df <- merge(df, df2[, .(patid, ckd_READ)])
df[, ckd_READ_bin := as.numeric(ckd_READ != 0)]
table(df$ckd_READ_bin, exclude = NULL)

#re-code smokers
table(df$smoking, exclude = NULL)
df[, smok_short := 0] # default value never-smokers; include those with NA
df[smoking == 2, smok_short := 1] # ex-smokers
df[smoking > 2, smok_short := 2] #  current smokers
table(df$smok_short, exclude = NULL)

#clinical CKD
df[, ckd_NICE := as.numeric(ckd != "ckd2")]
df[ckd == "ckd2" & (crd == 1 | acr > 1), ckd_NICE := 1]
table(df$ckd_NICE, exclude = NULL)

# update ckd_READ
df[, ckd_READ := as.numeric(as.character(ckd_READ))]
df[is.na(ckd_READ), ckd_READ := 0]
table(df$ckd_READ, exclude = NULL)

vars <- list(
  list(varName = "age", cts = TRUE, digits = 0),
  list(varName = "sex", cts = FALSE, digits = 0),
  list(varName = "smok_short", cts = FALSE, digits = 0),
  list(varName = "bmi", cts = TRUE, digits = 0),
  list(varName = "sbp", cts = TRUE, digits = 0),
  list(varName = "dbp", cts = TRUE, digits = 0),
  list(varName = "hdl", cts = TRUE, digits = 1),
  list(varName = "ldl", cts = TRUE, digits = 1),
  list(varName = "chol", cts = TRUE, digits = 1),
  list(varName = "hyp_diag", cts = FALSE, digits = 0),
  #list(varName = "cvd", cts = FALSE, digits = 0),
  list(varName = "diab_84", cts = FALSE, digits = 0),
  list(varName = "ra", cts=FALSE, digits=0),
  list(varName = "fh_chd", cts=FALSE, digits=0),
  list(varName = "cancer", cts = FALSE, digits = 0),
  list(varName = "ckd_NICE", cts = FALSE, digits = 0),
  list(varName = "ckd_READ_bin", cts = FALSE, digits = 0),
  list(varName = "acr", cts = FALSE, digits = 0),
  list(varName = "Base_Cardiovascular", cts=FALSE, digits=0),
  list(varName = "Base_Cerebrovascular", cts=FALSE, digits=0),
  list(varName = "Base_HF", cts=FALSE, digits=0),
  list(varName = "statin_84", cts = FALSE, digits = 0),
  list(varName = "antihyp", cts = FALSE, digits = 0),
  list(varName = "apc_base_84", cts = FALSE, digits = 0),
  list(varName = "imd", cts=FALSE, digits=0),
  list(varName = "region", cts=FALSE, digits=0)
)

# by CKD stage for all 1.13mm
alpha <- demographics(df = data.frame(df), vars = vars, subgroupVar = "ckd")
alpha
write.csv(alpha, file = paste(outputDir,  "demographics table (ALL).csv", sep = ""), row.names = FALSE)
table(df$ckd)

# by CKD stage for estimation cohort
alpha <- demographics(df = data.frame(df[estimation==1]), vars = vars, subgroupVar = "ckd")
alpha
write.csv(alpha, file = paste(outputDir,  "demographics table (estimation).csv", sep = ""), row.names = FALSE)
table(df$ckd)

# by CKD stage for validation cohort
alpha <- demographics(df = data.frame(df[estimation==0]), vars = vars, subgroupVar = "ckd")
alpha
write.csv(alpha, file = paste(outputDir,  "demographics table (validation).csv", sep = ""), row.names = FALSE)
table(df[estimation==0]$ckd)
table(df[estimation==1]$ckd)

#############################################################
### Endpoints
#############################################################

rm(list = ls())

library(tidyverse)
library(data.table)

source("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\src\\model_paper\\numbers_functions.R")

outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Documents\\Papers\\Model Paper\\outputs\\"

df <- data.table(get(load("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\data\\processed\\df_Aevents_cross_validation.Rdata")))

## Endpoints
vars2 <- list(
  list(varName = "NVD", cts = FALSE, digits = 0),
  list(varName = "VD", cts = FALSE, digits = 0),
  list(varName = "Stroke_or_VD", cts = FALSE, digits = 0),
  list(varName = "MI_Stroke_or_VD", cts = FALSE, digits = 0),
  list(varName = "HF_hosp", cts=FALSE, digits=0))

endpts <- df[,lapply(.SD, sum), .SDcols=c("NVD", "VD", "Stroke_or_VD", "MI_Stroke_or_VD", "HF_hosp"), by=c("patid", "ckd", "estimation")]
endpts

# everyone
beta <- demographics(df=data.frame(endpts), vars=vars2, subgroupVar = "ckd")
beta
write.csv(beta, file = paste(outputDir,  "Events table (All).csv", sep = ""), row.names = FALSE)

# estimation cohort
beta <- demographics(df=data.frame(endpts[estimation==1]), vars=vars2, subgroupVar = "ckd")
beta
write.csv(beta, file = paste(outputDir,  "Events table (Estimation).csv", sep = ""), row.names = FALSE)

# validation cohort
beta <- demographics(df=data.frame(endpts[estimation==0]), vars=vars2, subgroupVar = "ckd")
beta
write.csv(beta, file = paste(outputDir,  "Events table (Validation).csv", sep = ""), row.names = FALSE)


#############################################################
### Age distribution
#############################################################

rm(list = ls())
library(tidyverse)

df2 <- get(load("K:/SHARP_EE/Monitoring CKD CHF/Data/Derived CPRD data/df_incBase_hist_incACR.Rdata"))

test <- df2 %>%
  filter(age < 60 & ckd == "ckd2")

mean(test$age)
