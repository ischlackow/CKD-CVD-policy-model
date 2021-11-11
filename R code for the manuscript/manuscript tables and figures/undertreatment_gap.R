### Part 1: produce the file with baseline characteristics ----------------------------------------

### This program takes the population from df_baseline_cohort, and updates their baseline characteristics on 01/01/2014 ###
### The code is an adaptation of df_baseline_cohort.R, but only the variables required for treatment allocation are recorded
### treatment prescription is recorded for 01/01/2014-31/12/2014 rather than at baseline

######################################################################
### Setup directories and load files
######################################################################

rm(list = ls())

library(readstata13)
library(dplyr)
library(gdata)
library(data.table)
library(tidyverse)

# directories
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
sourceDir_R_raw <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"
sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_lookup <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\TXTFILES\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Documents\\Papers\\Model Paper\\aux data\\"

# functions
source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions v3.R")
source("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\src\\extrapolation\\preparation_functions.R" )
#source("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\src\\model_paper\\model_paper_functions.R" )

# original cohort
setwd(sourceDir_R)
df_base <- get(load("df_base_cohort.Rdata"))

# other datasets
setwd(sourceDir_R)
# information on date of birth
df_pat <- data.table(get(load("df_pat.Rdata")))
df_pat <- df_pat[, .(patid, dob)]
# eGFR tests
df_egfr <- get(load("df_egfr.Rdata"))

# CPRD lookups
setwd(sourceDir_lookup)
OPR <- data.table(read.delim("OPR.txt", stringsAsFactors = FALSE))
setkey(OPR, Code)
SUM <- data.table(read.delim("SUM.txt", stringsAsFactors = FALSE))
setkey(SUM, Code)
TQU <- data.table(read.delim("TQU.txt", stringsAsFactors = FALSE))
setkey(TQU, Code)
YND <- data.table(read.delim("YND.txt", stringsAsFactors = FALSE))
setkey(YND, Code)

######################################################################
### New baseline date
######################################################################

d_baseline <- as.Date("2013-01-01")

######################################################################
### select patients & record fixed variables
######################################################################

# keep only those in study on d_baseline
df_base <- df_base[studyentry <= d_baseline & fuenddate > d_baseline]
nrow(df_base)

range(df_base$studyentry)
range(df_base$fuenddate)
hist(as.numeric(df_base$fuenddate - df_base$studyentry))
hist(as.numeric(df_base$fuenddate - d_baseline))
quantile(as.numeric(df_base$fuenddate - d_baseline)) # 25% = 319, so the vas majority is almost a year

# keep fixed variables
df <- df_base[, .(patid, sex, ethnicity, imd, fh_chd, fuenddate)]
df[, studyentry := d_baseline]

# re-generate the imd variable
df[, imd_cts := NaN]
df[imd == 1, imd_cts := 4.51]
df[imd == 2, imd_cts := 0.84]
df[imd == 3, imd_cts := -1.05]
df[imd == 4, imd_cts := -2.17]
df[imd == 5, imd_cts := -3.15]
df[, imd := NULL]
table(df$imd_cts, exclude = NULL)

# define id_ckd to be used later
id_ckd <- df$patid

######################################################################
### update other baseline variables
######################################################################

### pick up the latest value up to d_baseline
### use the .latest function from function v2.R

### age ###

df_pat[, d := d_baseline]
df_pat[, age := round(as.numeric(d - dob) / 365.25)]
df <- merge(df, df_pat[, .(patid, age)], by = "patid", all.x = TRUE, sort = TRUE)

### CKD stage ###

#just need the latest eGFR measurement, as this is past patient's cohort entry
test <- .latest(df_egfr[eventdate < d_baseline])
df <- merge(df, test[, .(patid, ckd)], by = "patid", all.x = TRUE, sort = TRUE)

### crd ###

id <- unique(.clean_df(filename = "df_crd.Rdata")[eventdate <= d_baseline]$patid)
df[, crd := as.numeric(patid %in% id)]
table(df$crd, exclude = NULL)

### ckd_READ ###

# latest CKD 1
df_ckd1 <- .clean_df(filename = "df_ckd1.Rdata")[eventdate <= d_baseline]
df_ckd1[, stage := 1]
df_ckd1 <- df_ckd1[, .(patid, eventdate, stage)]
df_ckd1 <- .latest(df_ckd1)
# latest CKD 2
df_ckd2 <- .clean_df(filename = "df_ckd2.Rdata")[eventdate <= d_baseline]
df_ckd2[, stage := 2]
df_ckd2 <- df_ckd2[, .(patid, eventdate, stage)]
df_ckd2 <- .latest(df_ckd2)
# latest CKD 3A
df_ckd3A <- .clean_df(filename = "df_ckd3A.Rdata")[eventdate <= d_baseline]
df_ckd3A[, stage := "3A"]
df_ckd3A <- df_ckd3A[, .(patid, eventdate, stage)]
df_ckd3A <- .latest(df_ckd3A)
# CKD 3B
df_ckd3B <- .clean_df(filename = "df_ckd3B.Rdata")[eventdate <= d_baseline]
df_ckd3B[, stage := "3B"]
df_ckd3B <- df_ckd3B[, .(patid, eventdate, stage)]
df_ckd3B <- .latest(df_ckd3B)
# CKD 3
df_ckd3 <- .clean_df(filename = "df_ckd3.Rdata")[eventdate <= d_baseline]
df_ckd3[, stage := 3]
df_ckd3 <- df_ckd3[, .(patid, eventdate, stage)]
df_ckd3 <- .latest(df_ckd3)
# CKD 4
df_ckd4 <- .clean_df(filename = "df_ckd4.Rdata")[eventdate <= d_baseline]
df_ckd4[, stage := 4]
df_ckd4 <- df_ckd4[, .(patid, eventdate, stage)]
df_ckd4 <- .latest(df_ckd4)
# CKD 5
df_ckd5 <- .clean_df(filename = "df_ckd5.Rdata")[eventdate <= d_baseline]
df_ckd5[, stage := 5]
df_ckd5 <- df_ckd5[, .(patid, eventdate, stage)]
df_ckd5 <- .latest(df_ckd5)
# combine & find the latest diagnosis
df_ckd <- rbind(df_ckd1, df_ckd2, df_ckd3A, df_ckd3B, df_ckd3, df_ckd4, df_ckd5)
df_ckd <- .latest(df_ckd)
df_ckd[, ckd_READ := factor(stage, levels = c(1, 2, "3A", "3B", 3, 4, 5))]
# add to df
df <- merge(df, df_ckd[, .(patid, ckd_READ)], by = "patid", all.x = TRUE)
table(df$ckd_READ, df$ckd, exclude = NULL)
table(df_base$ckd_READ, df_base$ckd, exclude = NULL)

### ckd45_or_crd ###

df[, ckd45_or_crd := (ckd_READ %in% c("4", "5") | crd == 1)]

### ckd_NICE ###

df[, ckd_NICE := (!is.na(ckd_READ) | crd == 1)]
table(df$ckd_NICE, df$ckd_READ, exclude = NULL)
table(df$ckd_NICE, df$crd, exclude = NULL)

### acr ###

# calculate year since randomisation
df_studyentry <- df_base[, .(patid, studyentry)]
df_studyentry[, d := d_baseline]
df_studyentry[, y := ceiling(as.numeric((d - studyentry) / 365.25))]
df_studyentry <- df_studyentry[, .(patid, y)]

# numeric
df_alb <- .clean_df(sourceDir = sourceDir_R, filename = "df_alb_peryear_TESTNUM.Rdata")
df_alb <- merge(df_alb, df_studyentry, by = "patid", all.x = TRUE)
df_alb <- df_alb[eventyear <= y]
# first consider stage 2 or 3 only
alpha <- df_alb[stage >= 2]
# latest stage before study entry
alpha <- alpha[eventyear == max(eventyear), .SD, by = "patid"]
# add to df
alpha <- alpha[, .(patid, stage)]
colnames(alpha)[2] <- "album_TESTNUM"
df <- merge(df, alpha[, .(patid, album_TESTNUM)], all.x = TRUE, by = "patid")
# now differentiate between A1 & no measurement
beta <- df_alb[stage == 1]
id <- beta$patid
df[patid %in% id & is.na(album_TESTNUM), album_TESTNUM := 1]
# the rest would be with no measurements then
table(df$album_TESTNUM, exclude = NULL)
table(df$album_TESTNUM, df$ckd, exclude = NULL)

# categorical (including diagnoses)
df_alb <- .clean_df(sourceDir = sourceDir_R, filename = "df_alb_peryear_TESTCAT_CLIN.Rdata")
df_alb <- merge(df_alb, df_studyentry, by = "patid", all.x = TRUE)
df_alb <- df_alb[eventyear <= y]
# first consider stage 2 or 3 only
alpha <- df_alb[stage >= 2]
# latest stage before study entry
alpha <- alpha[eventyear == max(eventyear), .SD, by = "patid"]
# add to df
alpha <- alpha[, .(patid, stage)]
colnames(alpha)[2] <- "album_TESTCAT"
df <- merge(df, alpha[, .(patid, album_TESTCAT)], all.x = TRUE, by = "patid")
# now differentiate between A1 & no measurement
beta <- df_alb[stage == 1]
id <- beta$patid
df[patid %in% id & is.na(album_TESTCAT), album_TESTCAT := 1]
# the rest would be with no measurements then
table(df$album_TESTCAT, exclude = NULL)
table(df$album_TESTCAT, df$ckd, exclude = NULL)
table(df$album_TESTNUM, df$album_TESTCAT, exclude = NULL)

# combine both categories into one
df[, acr := album_TESTNUM]
df[is.na(acr) & !is.na(album_TESTCAT), acr := album_TESTCAT, ]
table(df$acr, exclude = NULL)
table(df$acr, df$album_TESTNUM, exclude = NULL)
table(df$acr, df$album_TESTCAT, exclude = NULL)
table(df$album_TESTCAT)
df[, album_TESTNUM := NULL]
df[, album_TESTCAT := NULL]
table(df$acr, exclude = NULL)
df[is.na(acr), acr := 0]
table(df$acr, exclude = NULL)

### smoking ###

df_smok <- .clean_df(filename = "df_smok.Rdata")[eventdate <= d_baseline]
df_smok <- .latest(df_smok)

# daily consumption
# convert data2
df_smok[, data2 := as.numeric(data2)]
df_smok[is.na(data2), data2 := 0]
df_smok[data2 >= 300, data2 := 0]
table(df_smok$data2, exclude = NULL)
# convert data3
df_smok[, data3 := as.numeric(data3)]
df_smok[is.na(data3), data3 := 0]
table(df_smok$data3, exclude = NULL)
# NB: Assume 1 cigar = 2 cigarettes. 
# Source: http://www.nhsggc.org.uk/about-us/professional-support-sites/cdm-local-enhanced-services/health-determinants/assess-status/smoking/
# convert data4
df_smok[, data4 := as.numeric(data4)]
df_smok[is.na(data4), data4 := 0]
table(df_smok$data4, exclude = NULL)
# NB Assume 1 oz of tobacco = 7 cigarettes per day
# Source: http://www.nhsggc.org.uk/about-us/professional-support-sites/cdm-local-enhanced-services/health-determinants/assess-status/smoking/
# add up
df_smok[, n_smok := data2 + 2 * data3 + 7 * data4]

# determine smoking status
# add codes
colnames(YND) <- c("data1", "status")
df_smok[, data1 := as.integer(data1)]
df_smok <- merge(df_smok, YND, by = "data1")
df_smok[, stopdate := as.Date(data6, format = "%d/%m/%Y")]
df_smok <- df_smok[, .(patid, n_smok, status, stopdate)]
df_smok <- df_smok[order(patid)]
# smoking category
df_smok[, category := NaN]
# never smokers
beta <- df_smok[n_smok == 0 & is.na(stopdate) & status %in% c("No", "Data Not Entered")]
id <- unique(beta$patid)
df_smok[patid %in% id, category := 1]
# ex-smokers
beta <- df_smok[n_smok == 0 & (!is.na(stopdate) | status == "Ex")]
# NB: this picks up those with stopdate & status = "Data Not Entered"
id <- unique(beta$patid)
df_smok[patid %in% id, category := 2]
# current smokers
beta <- df_smok[n_smok > 0 | (n_smok == 0 & is.na(stopdate) & status %in% c("Yes"))]
# light smokers
gamma <- beta[n_smok <= 9]
id <- unique(gamma$patid)
df_smok[patid %in% id, category := 3]
# moderate smokers
gamma <- beta[n_smok >= 10 & n_smok <= 19]
id <- unique(gamma$patid)
df_smok[patid %in% id, category := 4]
# heavy smokers
gamma <- beta[n_smok >= 20]
id <- unique(gamma$patid)
df_smok[patid %in% id, category := 5]
# check
table(df_smok$category, exclude = NULL)

# add to df
df_smok <- df_smok[, .(patid, category)]
colnames(df_smok)[2] <- "smoking"
df <- merge(df, df_smok, by = "patid", all.x = TRUE)
table(df$smoking, exclude = NULL)
table(df_base$smoking, exclude = NULL)

### bmi ###

# height dataset
df_height <- .clean_df(filename = "df_height.Rdata")
df_height <- df_height[medcode %in% c(3, 158, 8105, 12362, 15546, 26474, 30763, 35476, 41045, 47516)]
# keep useful information only
df_height[, cm := as.numeric(data1)]
df_height[, data1 := NULL]
df_height <- distinct(df_height)
# remove implausible values
df_height <- df_height[cm >= 0.5 & cm <= 2.30]
# take averages for readings taken on the same day
df_height <- .dailymean(df_height, "cm")
# extract entry closest to d_baseline
df_height[, d := d_baseline]
df_height[, diff := as.numeric(abs(eventdate - d))]
df_height <- df_height[order(patid, diff)]
df_height <- df_height[, head(.SD, 1), by = patid]
# take the earliest entry, in case there were two measurements, exactly the same number of days either side of study entry
df_height <- df_height[order(patid, eventdate)]
df_height <- df_height[, head(.SD, 1), by = patid]
df_height <- df_height[, .(patid, cm)]

# add to df
df <- merge(df, df_height, by = "patid", all.x = TRUE)
sum(is.na(df$cm)) / nrow(df)

# weight dataset
df_weight <- .clean_df(filename = "df_weight.Rdata")
df_weight <- df_weight[medcode %in% c(2, 126, 158, 430, 1018, 1581, 2839, 6545, 6713, 7984, 8041, 8105, 
                        13200, 13278, 16404, 21520, 22556, 23376, 25951, 26473, 28937, 28946, 29029, 29538, 
                        32914, 32974, 37937, 38632, 42309, 101043, 103499, 104002)]
# clean
df_weight[, c("kg", "bmi") := list(as.numeric(data1), as.numeric(data3))]
df_weight <- df_weight[, .(patid, eventdate, kg, bmi)]
df_weight <- df_weight[eventdate <= d_baseline]
# remove implausible values
df_weight[kg < 25, kg := NA]
df_weight[bmi < 10 | bmi > 300, bmi := NA]
# extract latest information
df_bmi <- df_weight[!is.na(bmi)]
df_bmi <- .latest(df_bmi)
df_bmi <- .dailymean(df_bmi, "bmi")
# weight
df_kg <- df_weight[!is.na(kg)]
df_kg <- .latest(df_kg)
df_kg <- .dailymean(df_kg, "kg")
# add height
df_cm <- df[, .(patid, cm)]
# extract the latest date of these two
df_diff <- merge(df_bmi, df_kg, by = "patid", all = TRUE, suffix = c("_bmi", "_kg"))
df_diff <- df_diff[, .(patid, eventdate_bmi, eventdate_kg)]
df_diff[, diff := as.numeric(eventdate_bmi - eventdate_kg)]

# add to df
# firstly, add those patients who have a bmi measurement, and no later weight measurement OR a bmi measurement only
alpha <- df_diff[eventdate_bmi >= eventdate_kg | is.na(eventdate_kg)]
id <- alpha$patid
beta <- df_bmi[patid %in% id]
df <- merge(df, beta[, .(patid, bmi)], by = "patid", all.x = TRUE)
sum(is.na(df$bmi)) / nrow(df)
# now, consider those with a bmi measurement and a later weight measurement OR a kg measurement ONLY
# so a later weight measurement will over-write a previous BMI measurement
alpha <- df_diff[eventdate_bmi < eventdate_kg | is.na(eventdate_bmi)]
id <- alpha$patid
beta <- df_kg[patid %in% id]
beta <- beta[, .(patid, kg)]
beta <- merge(beta, df_cm, by = "patid")
beta <- beta[!is.na(cm)]
beta[, bmi := kg / cm^2]
beta <- beta[bmi >= 10 & bmi <= 300]
id <- beta$patid
beta <- beta[, .(patid, bmi)]
colnames(beta)[2] <- "bmi_2"
df <- merge(df, beta, all.x = TRUE)
df[is.na(bmi), bmi := bmi_2]
df[, cm := NULL]
df[, bmi_2 := NULL]
sum(is.na(df$bmi)) / nrow(df)
quantile(df$bmi, na.rm = TRUE)
quantile(df_base$bmi, na.rm = TRUE)

### sbp ###

df_bp <- .clean_df(filename = "df_bp.Rdata")[eventdate <= d_baseline]
df_sbp <- df_bp[!(medcode %in% c(0, 2, 3, 6, 33, 34, 41, 47, 71, 76, 78, 89, 90, 93, 96, 
                                 103, 128, 134, 154, 161, 195, 241, 243, 252, 273, 337, 374, 412, 445, 451, 
                                 504, 532, 554, 630, 636, 637, 677, 693, 709, 711, 758, 805, 814, 844, 846, 984,
                                 1035, 1135, 1206, 1289, 1298, 1344, 1430, 1433, 1446, 1661, 1792, 1864, 1880, 
                                 2043, 2074, 2084, 2122, 2157, 2317, 2378, 2379, 2430, 2542, 2627, 3442, 3627, 
                                 5803, 5829, 6154, 6166, 6336, 6470, 6677, 6874, 7191, 7917, 8141, 9484,
                                 10043, 10558, 10632, 10674, 10711, 11427, 11717, 12350, 12566, 12936, 12947, 13144, 13185, 13228, 14704, 14805, 19506,
                                 20354, 20696, 22595, 22938, 22942, 28554,
                                 31305, 33330, 37243, 37312, 41445, 42280, 43547, 48008, 51357, 65990,
                                 100958))]
# keep useful information only
df_sbp[, sbp := as.numeric(data2)]
df_sbp <- df_sbp[, .(patid, eventdate, sbp)]
# remove duplicates
df_sbp <- distinct(df_sbp)
# leave plausible values only
df_sbp <- df_sbp[sbp >= 30 & sbp <= 350]
# take averages for readings taken on the same day
df_sbp <- .dailymean(df_sbp, "sbp")
# latest information
df_sbp_latest <- .latest(df_sbp)
# add to df
df <- merge(df, df_sbp_latest, by = "patid", all.x = TRUE)
quantile(df$sbp, na.rm = TRUE)
sum(is.na(df$sbp)) / nrow(df)

### dbp ###

df_dbp <- df_bp[!(medcode %in% c(0, 2, 3, 6, 33, 34, 41, 47, 71, 76, 78, 90, 93, 96, 
                                 128, 134, 154, 161, 241, 243, 252, 273, 337, 374, 445, 451, 
                                 504, 532, 554, 630, 636, 677, 693, 709, 711, 758, 805, 814, 844, 846, 984,
                                 1035, 1135, 1206, 1289, 1298, 1344, 1430, 1433, 1446, 1661, 1864, 1880, 
                                 2043, 2074, 2084, 2122, 2157, 2317, 2378, 2379, 2627, 3442, 3627, 
                                 5803, 5829, 6154, 6166, 6336, 6470, 6677, 6874, 7917, 8141, 9484,
                                 10043, 10558, 10632, 10674, 10711, 11427, 11717, 12350, 12947, 13144, 13185, 13228, 14704, 14805, 18418, 19506,
                                 20354, 20696, 22938, 22942, 28554, 29261,
                                 33330, 34618, 37242, 37312, 38277, 38278, 
                                 41052, 43282, 
                                 51357,
                                 100958))]

# keep useful information only
df_dbp[, dbp := as.numeric(data1)]
df_dbp <- df_dbp[, .(patid, eventdate, dbp)]
# leave plausible values only
df_dbp <- df_dbp[dbp >= 30 & dbp <= 350]
# remove duplicates
df_dbp <- distinct(df_dbp)
# take averages for readings taken on the same day
df_dbp <- .dailymean(df_dbp, colname = "dbp")
# latest information
df_dbp <- .latest(df_dbp)
# add to df
df <- merge(df, df_dbp[, .(patid, dbp)], by = "patid", all.x = TRUE)
quantile(df$dbp, na.rm = TRUE)
quantile(df_base$dbp, na.rm = TRUE)
sum(is.na(df$dbp)) / nrow(df)

### previous CV disease ###

# CPRD
df_cv <- .clean_df(filename = "df_cv.Rdata")
df_cv_pre <- df_cv[eventdate <= d_baseline]
df_cv_post <- df_cv[eventdate > d_baseline]
# HES
df_hes <- .clean_df(sourceDir = sourceDir_R, filename = "df_hes.Rdata", eventdate_colname = "discharged")
df_hes_pre <- df_hes[discharged <= d_baseline]
setwd(sourceDir_codelist)
codes_cv_ICD10 <- fread("CVD_ICD10_final.csv")
icd10_cv <- codes_cv_ICD10$ICD10_code
df_hes_pre[, ICD_3 := substr(ICD, 0, 3)]
df_hes_pre_cv <- df_hes_pre[ICD_3 %in% icd10_cv]
# combine
id <- unique(union(df_cv_pre$patid, df_hes_pre_cv$patid))
df[,  cvd := as.numeric(patid %in% id)]
table(df$cvd, exclude = NULL)
table(df_base$cvd, exclude = NULL)

### af ###

id <- unique(.clean_df(filename = "df_af.Rdata")[eventdate <= d_baseline]$patid)
df[, af := as.numeric(patid %in% id)]
table(df$af, exclude = NULL)
table(df_base$af, exclude = NULL)

### ra ###

id <- unique(.clean_df(filename = "df_ra.Rdata")[eventdate <= d_baseline]$patid)
df[, ra := as.numeric(patid %in% id)]
table(df$ra, exclude = NULL)
table(df_base$ra, exclude = NULL)

### tr_hyp_84 ###

# diagnosis
id_diag <- unique(.clean_df(filename = "df_hyp_diag.Rdata")[eventdate <= d_baseline]$patid)
# prescription of anti-hypertensives
id_Tx_84 <- .get_id_Tx_pred1(filename = "df_antihyp.Rdata", d1 = d_baseline, days = 84)
# combine
id <- intersect(id_diag, id_Tx_84)
df[, tr_hyp_84 := as.numeric(patid %in% id)]
table(df$tr_hyp_84, exclude = NULL)
table(df_base$tr_hyp_84, exclude = NULL)

### diab_84 ###

# diagnosis
# all pre-study entry records
df_diab <- .clean_df(filename = "df_diab_diag.Rdata")[eventdate <= d_baseline]
# extract first ever record; assuming this is the diagnosis date
df_diab[, medcode := NULL]
df_diab <- distinct(df_diab)
df_diab <- df_diab[order(patid, eventdate)]
df_diab <- df_diab[, head(.SD, 1), by = patid]
nrow(df_diab) == length(unique(df_diab$patid)) # yes, one record per patient
# calculate age at this date
df_diab <- merge(df_diab, df_pat[, .(patid, dob)], by = "patid")
df_diab[, age := as.numeric(df_diab$eventdate - df_diab$dob)/365.25]
# insulin treatment
id_Tx_84 <- .get_id_Tx_pred1(filename = "df_ins.Rdata", d1 = d_baseline, days = 84)
# combine information and determine type of diabetes
df_diab[, ins_84 := as.numeric(df_diab$patid %in% id_Tx_84)]
df_diab[, type_84 := ifelse(age <= 35 & df_diab$ins_84 == 1, 1, 2)]
# merge with df
df_diab <- df_diab[, .(patid, type_84)]
colnames(df_diab)[2] <- "diab_84"
df <- merge(df, df_diab, by = "patid", all.x = TRUE)
df[is.na(diab_84), diab_84 := 0]
table(df$diab_84, exclude = NULL)
table(df_base$diab_84, exclude = NULL)

## hdl cholesterol ###

df_hdl <- .clean_ctsvar(filename = "df_hdl.Rdata", colname = "hdl", pre = FALSE)[eventdate <= d_baseline]
df_hdl <- df_hdl[!(medcode %in% c(2379, 11867, 93756))]
df_hdl <- df_hdl[Specimen.Unit.Of.Measure %in% c("g/L", "mg/dL", "mmol/L", "mol/L")]
# convert into mmol/L
df_hdl[, hdl.mmol.L := NaN]
df_hdl[Specimen.Unit.Of.Measure == "g/L", hdl.mmol.L := hdl * 100 * 0.02586]
df_hdl[Specimen.Unit.Of.Measure == "mg/dL", hdl.mmol.L := hdl * 0.02586]
df_hdl[Specimen.Unit.Of.Measure == "mmol/L", hdl.mmol.L := hdl]
df_hdl[Specimen.Unit.Of.Measure == "mol/L", hdl.mmol.L := hdl * 1000]
df_hdl <- df_hdl[Specimen.Unit.Of.Measure != "mol/L"]
# remove implausible values
df_hdl <- df_hdl[hdl.mmol.L >= 0.1 & hdl.mmol.L <= 30]
# mean per day
df_hdl <- .dailymean(df_hdl, colname = "hdl.mmol.L")
# latest values
df_hdl <- .latest(df_hdl)

### total cholesterol ###

df_totchol <- .clean_ctsvar(filename = "df_serum_chol.Rdata", colname = "chol", pre = FALSE)[eventdate <= d_baseline]
df_totchol <- df_totchol[medcode %in% c(12, 622, 2493, 10940, 12821, 13733, 18040, 18147, 18443, 26902, 29202, 35720)]
df_totchol <- df_totchol[Specimen.Unit.Of.Measure %in% c("g/dL", "mg/dL", "mg/L", "mmol/L", "mol/L", "umol/L")]
df_totchol[, chol.mmol.L := NaN]
df_totchol[Specimen.Unit.Of.Measure == "g/dL", chol.mmol.L := chol * 1000 * 0.02586]
df_totchol[Specimen.Unit.Of.Measure == "mg/dL", chol.mmol.L := chol * 0.02586]
df_totchol[Specimen.Unit.Of.Measure == "mg/L", chol.mmol.L := chol * 0.1 * 0.02586]
df_totchol[Specimen.Unit.Of.Measure == "mmol/L", chol.mmol.L := chol]
df_totchol[Specimen.Unit.Of.Measure == "mol/L", chol.mmol.L := chol * 1000]
df_totchol <- df_totchol[Specimen.Unit.Of.Measure != "mol/L"]
df_totchol[Specimen.Unit.Of.Measure == "umol/L", chol.mmol.L := chol * 0.001]
df_totchol <- df_totchol[Specimen.Unit.Of.Measure != "umol/L"]
# remove implausible values
df_totchol <- df_totchol[chol.mmol.L >= 0.1 & chol.mmol.L <= 30]
# mean per day
df_totchol <- .dailymean(df_totchol, colname = "chol.mmol.L")
df_totchol <- .latest(df_totchol)

### chol_ratio ###

# Scenario 1: information on the ratio is available

# leave pre & post values
df_chol <- .clean_ctsvar(filename = "df_chol.Rdata", colname = "ratio", pre = FALSE)[eventdate <= d_baseline]
df_chol <- df_chol[medcode %in% c(14108, 14370, 14371, 14372, 40935)]
df_chol <- df_chol[Specimen.Unit.Of.Measure %in% c("1/1", "No Data Entered", "ratio")]
# convert into ratio; units are the same but numerators/denominators aren't
df_chol[, chol_hdl := NaN]
df_chol[medcode %in% c(14108), chol_hdl := 1 / ratio]
df_chol <- df_chol[medcode != 14108]
df_chol[medcode %in% c(14370), chol_hdl := (1 / ratio) + 1]
df_chol <- df_chol[medcode != 14370]
df_chol[medcode %in% c(14371, 14372, 40935), chol_hdl := ratio]
# filter
df_chol <- df_chol[chol_hdl >= 1 & chol_hdl <= 500] 
# clean & add to the output dataset
# take mean of tests taken on the same day
df_chol <- .dailymean(df_chol, colname = "chol_hdl")
df_ratio <- df_chol

# Scenario 2: There are separate tests for Total & HDL cholesterol

# combine df_hdl & df_totchol
df_chol <- merge(df_hdl, df_totchol, by = c("patid", "eventdate"))
df_chol[, chol_hdl := chol.mmol.L / hdl.mmol.L]
df_chol <- df_chol[, .(patid, eventdate, chol_hdl)]
# combine with df_ratio
df_chol <- rbind(df_ratio, df_chol)
df_chol <- distinct(df_chol)
# add dates of study entry, cv and end of follow up
# end of follow up incorporates cv date where relevant
df_chol <- merge(df_chol, df[, .(patid, fuenddate)], by = "patid")
df_chol <- df_chol[eventdate <= fuenddate]
# add date of statin prescription
# extract all precriptions
df_statin <- .clean_df(filename = "df_statins.Rdata")
# NB: at least two prescription would have been required at the extraction stage
# So just take the min date as the date of statin prescription
df_statin <- df_statin[order(patid, eventdate)]
df_statin <- df_statin[, head(.SD, 1), by = patid]
colnames(df_statin)[3] <- "statindate"
# add information to df_chol
df_chol <- merge(df_chol, df_statin[, .(patid, statindate)], by = "patid", all.x = TRUE)
df_chol[, censor := pmin(fuenddate, statindate, na.rm = TRUE)]
df_chol <- df_chol[eventdate <= censor]
df_chol <- df_chol[, .(patid, eventdate, chol_hdl)]
# extract value closest to  baseline
df_chol[, diff := as.numeric(abs(eventdate - d_baseline))]
df_chol <- df_chol[order(patid, diff)]
df_chol <- df_chol[, head(.SD, 1), by = patid]

# add to df
df <- merge(df, df_chol[, .(patid, chol_hdl)], by = "patid", all.x = TRUE)
sum(is.na(df$chol_hdl)) / nrow(df) # about 38%
colnames(df)[which(colnames(df) == "chol_hdl")] <- "chol_ratio"

quantile(df$chol_hdl, na.rm = TRUE)
quantile(df_base$chol_hdl, na.rm = TRUE)

### additional variables for QRISK3 #######

### lupus ###

id <- unique(.clean_df(filename = "df_lupus.Rdata")[eventdate <= d_baseline]$patid)
df[, lupus := as.numeric(patid %in% id)]
table(df$lupus, exclude = NULL)
table(df_base$lupus, exclude = NULL)

### antipsych_84

id <- .get_id_Tx_pred1(filename = "df_antipsych.Rdata", d1 = d_baseline, days = 84)
df[, antipsych_84 := as.numeric(patid %in% id)]
table(df$antipsych_84, exclude = NULL)
table(df_base$antipsych_84, exclude = NULL)

### steroids_84

id <- .get_id_Tx_pred1(filename = "df_steroids.Rdata", d1 = d_baseline, days = 84)
df[, steroids_84 := as.numeric(patid %in% id)]
table(df$steroids_84, exclude = NULL)
table(df_base$steroids_84, exclude = NULL)

### impotence_84

# diagnosis
id_diag <- unique(.clean_df(filename = "df_erectile_diag.Rdata")[eventdate <= d_baseline]$patid)

# treatment
id_Tx_84 <- .get_id_Tx_pred1(filename = "df_erectile_Tx.Rdata", d1 = d_baseline, days = 84)

# combine & add to df_base
id_84 <- union(id_diag, id_Tx_84)
df[, impotence_84 := as.numeric(patid %in% id_84)]
# remove for females
df[sex == "F", impotence_84 := 0]
table(df$impotence_84, exclude = NULL)

### migraine

id <- unique(.clean_df(filename = "df_migraine.Rdata")[eventdate <= d_baseline]$patid)
df[, migraine := as.numeric(patid %in% id)]
table(df$migraine, exclude = NULL)
table(df_base$migraine, exclude = NULL)

### mental

id <- unique(.clean_df(filename = "df_mental.Rdata")[eventdate <= d_baseline]$patid)
df[, mental := as.numeric(patid %in% id)]
table(df$mental, exclude = NULL)
table(df_base$mental, exclude = NULL)

### sbp_sd

# all BP values in the five years before the baseline date
df_t <- merge(df_sbp, df_studyentry, by = "patid")
df_t <- df_t[d_baseline - eventdate <= 365.25 * 5]
df_t <- distinct(df_t)
nrow(df_t)

# leave only patients with at least two measurements
alpha <- df_t[, .N, by = patid]
alpha <- alpha[ N >= 2]
id <- alpha$patid
df_t <- df_t[patid %in% id]

# calculate sd
df_t <- df_t[, .(sbp_sd = sd(sbp)), by = patid]

# add to df_base
df <- merge(df, df_t, by = "patid", all.x = TRUE)
quantile(df$sbp_sd, na.rm = TRUE)
# replace 0 with 0.1 [checks that this seems okay performed earlier, during Pengfei's thesis]
# this does indeed seem to correspond to patients whose records all contained the same value
# this could be genuine, a misprint, or perhaps a value carried over from the previous visit
# but we cannot tell so assume correct as is.
df[sbp_sd == 0, sbp_sd := 0.1]
quantile(df$sbp_sd, na.rm = TRUE)
quantile(df_base$sbp_sd, na.rm = TRUE)
sum(is.na(df$sbp_sd)) / nrow(df) # 7% missing
nrow(df) == length(unique(df$patid))

# save
setwd("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\data\\processed\\")
save(df, file = "df_2013_pre_tx.Rdata")

### Part 2: Eligibility for treatment -----------------------------------

df[, ckd_READ := as.character(ckd_READ)]
df[is.na(ckd_READ), ckd_READ := "none"]
table(df$ckd_READ, exclude = NULL)

### calculate QRISK3 score

.temp_clean_df_qrisk <- function(df_t){
  df_t <- df_t %>%
    # use imputed variables while over-writing old variables if required
    rename(town = imd_cts, rati = chol_ratio) %>%
    rename(secondary = cvd) %>%
    # other
    rename(b_sle = lupus, b_treatedhyp = tr_hyp_84,
           b_AF =  af, b_atypicalantipsy = antipsych_84, b_corticosteroids = steroids_84, b_impotence2 = impotence_84,
           b_migraine = migraine, b_ra = ra, b_semi = mental) %>% 
    mutate(b_type1 = as.numeric(diab_84 == 1), b_type2 = as.numeric(diab_84 == 2)) %>%
    rename(fh_cvd = fh_chd, smoke_cat = smoking) %>%
    #mutate(sbp = sbp + 139) %>%
    #mutate(b_renal = as.numeric(acr > 1 | crd == 1 | (ckd != "ckd2" & ckd != "ckd1"))) %>%
    mutate(b_renal = (ckd_READ %in% c("3", "3A", "3B", "4", "5"))) %>%
    mutate(ethnicity = ifelse(ethnicity == 11 | is.na(ethnicity), 1, ethnicity)) %>%
    mutate(sbps5 = ifelse(is.na(sbp_sd), mean(sbp_sd, na.rm=TRUE), sbp_sd)) %>%
    select(c(patid, secondary, age, sex, ckd_READ, ethnicity, town,
             smoke_cat, b_type1, b_type2, fh_cvd, b_renal, b_AF, b_treatedhyp, b_migraine, b_ra, b_sle, b_semi, b_atypicalantipsy, b_corticosteroids,
             b_impotence2, rati, sbp, sbps5, bmi, acr, dbp)) %>%
    arrange(patid)
  
  return(df_t)
  
}

# females
alpha <- .temp_clean_df_qrisk(filter(df, sex == "F"))
alpha <- mutate(alpha, sbp_orig = sbp)
df_female <- Qrisk3_female(alpha)
hist(df_female$score)

# males
beta <- .temp_clean_df_qrisk(filter(df, sex == "M"))
beta <- mutate(beta, sbp_orig = sbp)
df_male <- Qrisk3_male(beta)
hist(df_male$score)

# combine
df_tx <- data.table(rbind(df_male, df_female))

# any ckd flag
df_tx$ckd_read_bin <- df_tx$ckd_READ != "none" 

threshold_statin <- 10
threshold_antihyp <- 10

# statins
df_tx[, need_statins := 0]
# Qrisk3 threshold
df_tx[score >= threshold_statin, need_statins := 1]
# Clinical CKD
df_tx[ckd_read_bin == 1, need_statins := 1]
# >= 85 years old
df_tx[age >= 85, need_statins := 1]
# type 1 diabetes
df_tx[b_type1 == 1, need_statins := 1]
# secondary prevention
df_tx[secondary == 1, need_statins := 1]

# antihypertensives
df_tx[, need_antihyp := 0]
# CKD and [(diabetes and micro/macroalbuminuria) or macroalbuminuria]
df_tx[ckd_read_bin == 1 & (b_type1 == 1 | b_type2 == 1) & acr > 1, need_antihyp := 1]
df_tx[ckd_read_bin == 1 & acr == 3, need_antihyp := 1]
# stage 1 and (CVD, CKD, diabetes or Qrisk > 10)
df_tx[(sbp_orig >= 140 | dbp >= 90) & secondary == 1, need_antihyp := 1]
df_tx[(sbp_orig >= 140 | dbp >= 90) & ckd_read_bin == 1, need_antihyp := 1]
df_tx[(sbp_orig >= 140 | dbp >= 90) & (b_type1 == 1 | b_type2 == 1), need_antihyp := 1]
df_tx[(sbp_orig >= 140 | dbp >= 90) & score >= threshold_antihyp, need_antihyp := 1]
# stage 2
df_tx[sbp_orig >= 160 | dbp >= 100, need_antihyp := 1]

quantile(df_tx$score, na.rm = TRUE)
quantile(df_tx$sbp_orig, na.rm = TRUE)
quantile(df_tx$dbp, na.rm = TRUE)

# antiplatelets
df_tx[, need_antip := 0]
df_tx[secondary == 1, need_antip := 1]

# Part 3: Treatment prescription -------------------------------------

d1 <- d_baseline
d2 <- d_baseline + 365.25
days <- 84

# statins
id <- .get_id_Tx_d1d2(filename = "df_statins.Rdata", d1 = d1, d2 = d2, days = days)
df_tx[, statin_84 := as.numeric(patid %in% id)]
table(df_tx$statin_84, exclude = NULL)
table(df_tx$statin_84, df_tx$need_statins, exclude = NULL)

# antihyp
id <- .get_id_Tx_d1d2(filename = "df_antihyp.Rdata", d1 = d1, d2 = d2, days = days)
df_tx[, antihyp_84 := as.numeric(patid %in% id)]
table(df_tx$antihyp_84, exclude = NULL)
table(df_tx$antihyp_84, df_tx$need_antihyp, exclude = NULL)

# apc
id <- .get_id_Tx_d1d2(filename = "df_apc.Rdata", d1 = d1, d2 = d2, days = days)
df_tx[, apc_84 := as.numeric(patid %in% id)]
table(df_tx$apc_84, exclude = NULL)
table(df_tx$apc_84, df_tx$need_antip, exclude = NULL)

### Part 4: Clean and save ----------------------------------

setwd("K:\\SHARP_EE\\Monitoring CKD CHF\\ckd-monitoring\\data\\processed\\")
save(df_tx, file = "df_tx_2013.Rdata")

### Part 5: calculate undertreatment gap -----------------------

temp <- filter(data.frame(df_tx), (secondary == 0 & !is.na(score)) | secondary == 1)
round(1 - nrow(filter(df_tx, is.na(score) & secondary == 0)) / nrow(filter(df_tx, secondary == 0)), digits = 2)

### statins

alpha <- select(temp, secondary, score, need_statins, statin_84)
alpha <- filter(alpha, need_statins == 1)

# primary
alpha0 <- filter(alpha, secondary == 0)
nrow(alpha0)
round(table(alpha0$statin_84, exclude = NULL) / nrow(alpha0), digits = 2)

# secondary
alpha1 <- filter(alpha, secondary == 1)
nrow(alpha1)
round(table(alpha1$statin_84, exclude = NULL) / nrow(alpha1), digits = 2)

### antihypertensives

alpha <- select(temp, secondary, score, need_antihyp, antihyp_84)
alpha <- filter(alpha, need_antihyp == 1)

# primary
alpha0 <- filter(alpha, secondary == 0)
nrow(alpha0)
round(table(alpha0$antihyp_84, exclude = NULL) / nrow(alpha0), digits = 2)

# secondary
alpha1 <- filter(alpha, secondary == 1)
nrow(alpha1)
round(table(alpha1$antihyp_84, exclude = NULL) / nrow(alpha1), digits = 2)

### antiplatelets

alpha <- select(temp, secondary, score, need_antip, apc_84)
alpha <- filter(alpha, need_antip == 1)

# primary
alpha0 <- filter(alpha, secondary == 0)
nrow(alpha0)
round(table(alpha0$apc_84, exclude = NULL) / nrow(alpha0), digits = 2)

# secondary
alpha1 <- filter(alpha, secondary == 1)
nrow(alpha1)
round(table(alpha1$apc_84, exclude = NULL) / nrow(alpha1), digits = 2)

### mean proportions

mean(c(29, 44))

mean(c(73, 64, 76))
