######################################################################
### Patients with at least two creatinine tests
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)
library(data.table)
library(dtplyr)

sourceDir_codelist <- "L:\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "D:\\CKD\\source\\"
sourceDir_R <- "L:\\Monitoring CKD CHF\\Data\\derived CPRD data\\"
outputDir <- "L:\\Monitoring CKD CHF\\Data\\derived CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### load the lookups

setwd(sourceDir_R)
entity <- get(load("entity.Rdata"))

### medcodes ###

setwd(sourceDir_codelist)
codes_t <- read.csv("serum_creatinine_final.csv")
medcodes_t <- codes_t$medcode[codes_t$included == 1]

# entity
grep("creat", entity$description, ignore.case = TRUE, value = TRUE)
subset(entity, description == "Serum creatinine")
enttype_t <- 165 # as per the ISAC protocol

### extract data

df_output <- NULL

setwd(paste(sourceDir_CPRD, "Test", sep = "\\"))

startdate <- as.Date("2004-07-01")

i <- 0
p <- link_pracid[1]

for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
 
  # load the practice dataset
  f <- paste("p", p, "_Test_cut.dta", sep = "")
  df_test <- read.dta13(f)
  df_test <- select(df_test, patid, eventdate, medcode, enttype, data1, data2, data3)
  df_test <- mutate(df_test, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  df_test <- setDT(df_test)
  
  ### initialise output and do basic cleaning
  
  alpha <- df_test
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove tests with empty value
  alpha <- filter(alpha, !is.na(data2))
  # remove duplicates
  alpha <- distinct(alpha)
  
  ### filter out creatinine tests
  
  # leave only relevant SCr values
  alpha <- filter(alpha, medcode %in% medcodes_t & enttype == enttype_t) 
  alpha <- select(alpha, patid, eventdate, medcode, data1, data2, data3)
  # NB: PC already only left patients who have been registered at least 1 year after uts, so no need to do further filtering on the left
  # ie even tests on the dat of crd are okay
  
  ### leave only two latest Creatinine tests >= 90 days apart before startdate (& all tests in between)
  
  # sort by patid & date
  alpha <- alpha[order(alpha$patid, alpha$eventdate), ] 
  # leave all tests after startdate
  alpha_post <- filter(alpha, eventdate >= startdate)
  # leave only two latest Creatinine tests >= 90 days apart before startdate
  alpha_pre <- filter(alpha, eventdate < startdate)
  
  if (nrow(alpha_pre) > 0) {
    # latest test
    beta_pre <- group_by(alpha_pre, patid)
    beta_pre <- summarise(beta_pre, eventdate = max(eventdate))
    # first test at least 90 days before that
    alpha_pre <- alpha_pre[, tdiff := eventdate - max(eventdate), by = patid]
    beta_pre <-  group_by(alpha_pre, patid)
    beta_pre <- filter(beta_pre, tdiff <= -90)
    beta_pre <- summarise(beta_pre, tdiff_max = max(tdiff))
    # keep this and all tests in between
    alpha_pre <- left_join(alpha_pre, beta_pre)
    alpha_pre <- filter(alpha_pre, tdiff >= tdiff_max)
  }
  # combine
  alpha_post <- ungroup(alpha_post)
  alpha <- rbind(alpha_post, alpha_pre)
  alpha <- alpha[order(alpha$patid, alpha$eventdate), ] 
  
  ### leave patients with at least two SCr measurements
  
  alpha <- group_by(alpha, patid)
  freq <- summarise(alpha, n = n())
  ids <- freq$patid[freq$n >= 2]
  alpha <- ungroup(alpha)
  alpha <- filter(alpha, patid %in% ids)
  print(length(unique(alpha$patid)))

  ### add to the output
  
  alpha <- mutate(alpha, pracid = p)
  df_output[[p]] <- alpha
  
}

print(i)
# combine
df_Scr <- do.call(rbind , df_output)

setwd(outputDir)
save(df_scr, file = "df_scr.Rdata")

######################################################################
### Albuminuria
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### eligible patients ###
#
#setwd(sourceDir_R)
#df_basic <- get(load("df_basic.Rdata"))
#df_basic <- select(df_basic, patid, start, end)

### medcodes ###

setwd(sourceDir_codelist)

## ACR
#codes_t <- read.csv("acr_final.csv")
#medcodes_acr <- codes_t$medcode[codes_t$included == 1]
#df_output_acr <- NULL

## PCR
#codes_t <- read.csv("pcr_final.csv")
#medcodes_pcr <- codes_t$medcode[codes_t$included == 1]
#df_output_pcr <- NULL

## Albumin
#codes_t <- read.csv("albumin_final.csv")
#medcodes_alb <- codes_t$medcode[codes_t$included == 1]
#df_output_alb <- NULL

# Protein 
codes_t <- read.csv("protein_final.csv")
medcodes_prot <- codes_t$medcode[codes_t$included == 1]
df_output_prot <- NULL

## Creatinine
#codes_t <- read.csv("urine_creatinine_final.csv")
#medcodes_cr <- codes_t$medcode[codes_t$included == 1]
#df_output_cr <- NULL

##AER
#codes_t <- read.csv("aer_final.csv")
#medcodes_aer <- codes_t$medcode[codes_t$included == 1]
#df_output_aer <- NULL

##PER
#codes_t <- read.csv("per_final.csv")
#medcodes_per <- codes_t$medcode[codes_t$included == 1]
#df_output_per <- NULL

### extract the data ###

#p <- link_pracid[1]

setwd(paste(sourceDir_CPRD, "Test", sep = "\\"))
i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  # load the practice dataset
  f <- paste("p", p, "_Test_cut.dta", sep = "")
  df_test <- read.dta13(f)
  df_test <- select(df_test, patid, eventdate, medcode, enttype, data1, data2, data3)
  df_test <- mutate(df_test, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  ### initialise output and do basic cleaning
  
  alpha <- df_test
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove duplicates
  alpha <- distinct(alpha)
  ## filter by patients & censor dates

  ## ACR
  #df <- filter(alpha, medcode %in% medcodes_acr)
  #df_output_acr[[p]] <- df
  #print(paste("ACR:", nrow(df)))
  
  ## PCR
  #df <- filter(alpha, medcode %in% medcodes_pcr)
  #df_output_pcr[[p]] <- df
  #print(paste("PCR:", nrow(df)))
  
  ## Albumin
  #df <- filter(alpha, medcode %in% medcodes_alb)
  #df_output_alb[[p]] <- df
  #print(paste("Albumin:", nrow(df)))
  
  # Protein
  df <- filter(alpha, medcode %in% medcodes_prot)
  df_output_prot[[p]] <- df
  print(paste("Protein:", nrow(df)))
  
  ## Creatinine
  #df <- filter(alpha, medcode %in% medcodes_cr)
  #df_output_cr[[p]] <- df
  #print(paste("Creatinine:", nrow(df)))
  
  ## AER
  #df <- filter(alpha, medcode %in% medcodes_aer)
  #df_output_aer[[p]] <- df
  #print(paste("AER:", nrow(df)))
  
  ## PER
  #df <- filter(alpha, medcode %in% medcodes_per)
  #df_output_per[[p]] <- df
  #print(paste("PER:", nrow(df)))
}

## combine
#df_acr <- do.call(rbind, df_output_acr)
#print(nrow(df_acr))
#df_pcr <- do.call(rbind, df_output_pcr)
#print(nrow(df_pcr))
#df_alb <- do.call(rbind, df_output_alb)
#print(nrow(df_alb))
df_prot <- do.call(rbind, df_output_prot)
print(nrow(df_prot))
#df_cr <- do.call(rbind, df_output_cr)
#print(nrow(df_cr))
#df_aer <- do.call(rbind, df_output_aer)
#print(nrow(df_aer))
#df_per <- do.call(rbind, df_output_per)
#print(nrow(df_per))

setwd(outputDir)
#save(df_acr, file = "df_acr.Rdata")
#save(df_pcr, file = "df_pcr.Rdata")
#save(df_alb, file = "df_alb.Rdata")
save(df_prot, file = "df_prot.Rdata")
#save(df_cr, file = "df_cr.Rdata")
#save(df_aer, file = "df_aer.Rdata")
#save(df_per, file = "df_per.Rdata")

######################################################################
### QRISK variables: Clinical
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### medcodes ###

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")
colname <- "medcode"

## CV disease
#medcodes_cv <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "CVD_final.csv", colname = colname)
#df_output_cv <- NULL

## ethnicity (CPRD)
#medcodes_ethnicity <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "ethnicity_final.csv", colname = colname)
#df_output_ethnicity <- NULL

## FH CHD
#medcodes_fh_chd <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "FH_CHD_final.csv", colname = colname)
#df_output_fh_chd <- NULL

## diabetes
#medcodes_diab_diag <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "diabetes_final.csv", colname = colname)
#df_output_diab_diag <- NULL

## hypertension
#medcodes_hyp_diag <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "hypertension_final.csv", colname = colname)
#df_output_hyp_diag <- NULL

## RA
#medcodes_ra <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "RA_final.csv", colname = colname)
#df_output_ra <- NULL

## CKD 4 or 5
#medcodes_ckd45 <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "CKD45_final.csv", colname = colname)
#df_output_ckd45 <- NULL

# CKD 3
medcodes_ckd3 <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "CKD3_final.csv", colname = colname)
df_output_ckd3 <- NULL

## CRD
#medcodes_crd <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "crd_final.csv", colname = colname)
#df_output_crd <- NULL

## migraine
#medcodes_migraine <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "migraine_final.csv", colname = colname)
#df_output_migraine <- NULL

## lupus
#medcodes_lupus <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "lupus_final.csv", colname = colname)
#df_output_lupus <- NULL

## mental illness
#medcodes_mental <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "mental_final.csv", colname = colname)
#df_output_mental <- NULL

## HIV or AIDS
#medcodes_hiv <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "hiv_final.csv", colname = colname)
#df_output_hiv <- NULL

## erectile dysfunction
#medcodes_erectile_diag <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "erectile_diag_final.csv", colname = colname)
#df_output_erectile_diag <- NULL

## mental illness
#medcodes_mental <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "mental_final.csv", colname = colname)
#df_output_mental <- NULL

## FH CKD
#medcodes_fh_ckd <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "FH_CKD_final.csv", colname = colname)
#df_output_fh_ckd <- NULL

## Kidney stones
#medcodes_kidney_stones <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "kidney_stones_final.csv", colname = colname)
#df_output_kidney_stone <- NULL

## Atrial Fibrillation
#medcodes_af <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "AF_final.csv", colname = colname)
#df_output_af <- NULL

### extract data

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  # load the practice dataset: clinical
  setwd(paste(sourceDir_CPRD, "Clinical", sep = "\\"))
  f <- paste("p", p, "_Clinical_cut.dta", sep = "")
  df_clinical <- read.dta13(f)
  df_clinical <- select(df_clinical, patid, eventdate, medcode)
  df_clinical <- mutate(df_clinical, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  # load the practice dataset: referral
  setwd(paste(sourceDir_CPRD, "Referral", sep = "\\"))
  f <- paste("p", p, "_Referral_cut.dta", sep = "")
  df_ref <- read.dta13(f)
  df_ref <- select(df_ref, patid, eventdate, medcode)
  df_ref <- mutate(df_ref, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  ### initialise output and do basic cleaning
  
  alpha <- rbind(df_clinical, df_ref)
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove duplicates
  alpha <- distinct(alpha)
  
  ## CV
  #df <- filter(alpha, medcode %in% medcodes_cv)
  #df_output_cv[[p]] <- df
  #print(paste("CV:", nrow(df)))
  
  ## ethnicity
  #df <- filter(alpha, medcode %in%  medcodes_ethnicity)
  #df_output_ethnicity[[p]] <- df
  #print(paste("ethnicity:", nrow(df)))
  
  ## FH CHD
  #df <- filter(alpha, medcode %in%  medcodes_fh_chd)
  #df_output_fh_chd[[p]] <- df
  #print(paste("FH CHD:", nrow(df)))
  
  ## diabetes
  #df <- filter(alpha, medcode %in%  medcodes_diab_diag)
  #df_output_diab_diag[[p]] <- df
  #print(paste("Diabetes:", nrow(df)))
  
  ## hypertension
  #df <- filter(alpha, medcode %in%  medcodes_hyp_diag)
  #df_output_hyp_diag[[p]] <- df
  #print(paste("Hypertension:", nrow(df)))
  
  ## RA
  #df <- filter(alpha, medcode %in%  medcodes_ra)
  #df_output_ra[[p]] <- df
  #print(paste("RA:", nrow(df)))
  
  ## CKD 4 or 5
  #df <- filter(alpha, medcode %in%  medcodes_ckd45)
  #df_output_ckd45[[p]] <- df
  #print(paste("CKD 4 or 5:", nrow(df)))
  
  # CKD 3
  df <- filter(alpha, medcode %in%  medcodes_ckd3)
  df_output_ckd3[[p]] <- df
  print(paste("CKD 3:", nrow(df)))
  
  ## CRD
  #df <- filter(alpha, medcode %in%  medcodes_crd)
  #df_output_crd[[p]] <- df
  #print(paste("CRD:", nrow(df)))
  
  ## migraine
  #df <- filter(alpha, medcode %in%  medcodes_migraine)
  #df_output_migraine[[p]] <- df
  #print(paste("Migraine:", nrow(df)))
  
  ## lupus
  #df <- filter(alpha, medcode %in%  medcodes_lupus)
  #df_output_lupus[[p]] <- df
  #print(paste("Lupus:", nrow(df)))
  
  ## mental illness
  #df <- filter(alpha, medcode %in%  medcodes_mental)
  #df_output_mental[[p]] <- df
  #print(paste("Mental illness:", nrow(df)))
  
  ## HIV or AIDS
  #df <- filter(alpha, medcode %in%  medcodes_hiv)
  #df_output_hiv[[p]] <- df
  #print(paste("HIV or AIDS:", nrow(df)))
  
  ## erectile dysfunction
  #df <- filter(alpha, medcode %in%  medcodes_erectile_diag)
  #df_output_erectile_diag[[p]] <- df
  #print(paste("Erectile dysfunction:", nrow(df)))
  
  ## mental illness
  #df <- filter(alpha, medcode %in%  medcodes_mental)
  #df_output_mental[[p]] <- df
  
  ## FH_CKD
  #df <- filter(alpha, medcode %in%  medcodes_fh_ckd)
  #df_output_fh_ckd[[p]] <- df
  
  ## kidney stones
  #df <- filter(alpha, medcode %in%  medcodes_kidney_stones)
  #df_output_kidney_stone[[p]] <- df
  
  ## AF
  #df <- filter(alpha, medcode %in% medcodes_af)
  #df_output_af[[p]] <- df
  #print(paste("AF:", nrow(df)))

}

# combine and save

setwd(outputDir)

#df_cv <- do.call(rbind, df_output_cv)
#print(nrow(df_cv))
#save(df_cv, file = "df_cv.Rdata")

#df_ethnicity <- do.call(rbind, df_output_ethnicity)
#print(nrow(df_ethnicity))
#save(df_ethnicity, file = "df_ethnicity.Rdata")

#df_fh_chd <- do.call(rbind, df_output_fh_chd)
#print(nrow(df_fh_chd))
#save(df_fh_chd, file = "df_fh_chd.Rdata")

#df_diab_diag <- do.call(rbind, df_output_diab_diag)
#print(nrow(df_diab_diag))
#save(df_diab_diag, file = "df_diab_diag.Rdata")

#df_hyp_diag <- do.call(rbind, df_output_hyp_diag)
#print(nrow(df_hyp_diag))
#save(df_hyp_diag, file = "df_hyp_diag.Rdata")

#df_ra <- do.call(rbind, df_output_ra)
#print(nrow(df_ra))
#save(df_ra, file = "df_ra.Rdata")

#df_ckd45 <- do.call(rbind, df_output_ckd45)
#print(nrow(df_ckd45))
#save(df_ckd45, file = "df_ckd45_all.Rdata")

df_ckd3 <- do.call(rbind, df_output_ckd3)
print(nrow(df_ckd3))
save(df_ckd3, file = "df_ckd3_all.Rdata")

#df_crd <- do.call(rbind, df_output_crd)
#print(nrow(df_crd))
#save(df_crd, file = "df_crd.Rdata")

#df_migraine <- do.call(rbind, df_output_migraine)
#print(nrow(df_migraine))
#save(df_migraine, file = "df_migraine.Rdata")

#df_lupus <- do.call(rbind, df_output_lupus)
#print(nrow(df_lupus))
#save(df_lupus, file = "df_lupus.Rdata")

#df_mental <- do.call(rbind, df_output_mental)
#print(nrow(df_mental))
#save(df_mental, file = "df_mental.Rdata")

#df_hiv <- do.call(rbind, df_output_hiv)
#print(nrow(df_hiv))
#save(df_erectile_diag, file = "df_erectile_diag.Rdata")

#df_erectile_diag <- do.call(rbind, df_output_erectile_diag)
#print(nrow(df_erectile_diag))
#save(df_hiv, file = "df_hiv.Rdata")

#df_mental <- do.call(rbind, df_output_mental)
#print(nrow(df_mental))
#save(df_mental, file = "df_mental.Rdata")

#df_fh_ckd <- do.call(rbind, df_output_fh_ckd)
#print(nrow(df_fh_ckd))
#save(df_fh_ckd, file = "df_fh_ckd.Rdata")

#df_fh_kidney_stones <- do.call(rbind, df_output_kidney_stones)
#print(nrow(df_kidney_stones))
#save(df_kidney_stones, file = "df_kidney_stones.Rdata")

#df_af <- do.call(rbind, df_output_af)
#print(nrow(df_af))
#save(df_af, file = "df_af.Rdata")

######################################################################
### QRISK variables: Test file
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)

sourceDir_codelist <- "L:\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "D:\\CKD\\source\\"
sourceDir_R <- "L:\\Monitoring CKD CHF\\Data\\derived CPRD data\\"
outputDir <- "L:\\Monitoring CKD CHF\\Data\\derived CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### entity

setwd(sourceDir_R)
load("entity.Rdata")

### entity types

# Chol/HDL ratio
enttype_chol <- subset(entity, description == "HDL/LDL ratio")$enttype # 338
df_output_chol <- NULL

# Serum cholesterol
enttype_chol <- subset(entity, description == "Serum cholesterol")$enttype # 163
df_output_chol <- NULL

# HDL cholesterol
enttype_hdl <- subset(entity, description == "High density lipoprotein")$enttype # 175
df_output_hdl <- NULL

# LDL cholesterol
enttype_ldl <- subset(entity, description == "Low density lipoprotein")$enttype # 177
df_output_ldl <- NULL

# Triglycerides
enttype_trig <- subset(entity, description == "Triglycerides")$enttype # 202
df_output_trig <- NULL

# Hba1c
enttype_hba1c <- subset(entity, description == "HbA1c - diabetic control")$enttype # 275
df_output_hba1c <- NULL

# Triglycerides
enttype_trig <- subset(entity, description == "Triglycerides")$enttype # 202
df_output_trig <- NULL

# Haemoglobin
enttype_haem <- subset(entity, description == "Haemoglobin")$enttype # 173
df_output_haem <- NULL

# Calcium
enttype_calc <- subset(entity, description == "Calcium adjusted")$enttype # 159
df_output_calc <- NULL

# ALT
enttype_alt <- subset(entity, description == "Alanine aminotransferase")$enttype # 155
df_output_alt <- NULL

# ALP
enttype_alp <- subset(entity, description == "Alkaline phosphatase")$enttype # 153
df_output_alp <- NULL

# Bilirubin
enttype_bil <- subset(entity, description == "Bilirubin")$enttype # 158
df_output_bil <- NULL

# Phosphate
enttype_phos <- subset(entity, description == "Serum inorganic phosphate")$enttype # 188
df_output_phos <- NULL

# Albumin
enttype_alb <- subset(entity, description == "Albumin")$enttype # 152
df_output_alb <- NULL

### extract data
i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  ### Test file ###
  
  setwd(paste(sourceDir_CPRD, "Test", sep = "\\"))
  
  # load the practice dataset: test
  f <- paste("p", p, "_Test_cut.dta", sep = "")
  df_test <- read.dta13(f)
  df_test <- select(df_test, patid, eventdate, medcode, enttype, data1, data2, data3)
  df_test <- mutate(df_test, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  # initialise output and do basic cleaning
  alpha <- df_test
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  ## remove tests with empty value <- NOT HERE, as tests may have different formats!
  #alpha <- filter(alpha, !is.na(data2))
  # remove duplicates
  alpha <- distinct(alpha)
  
  # Cholesterol
  df <- filter(alpha, enttype %in% enttype_chol)
  df_output_chol[[p]] <- df
  print(paste("Chol:", nrow(df)))
  
  # Cholesterol
  df <- filter(alpha, enttype %in% enttype_chol)
  df_output_chol[[p]] <- df
  
  # HDL cholesterol
  df <- filter(alpha, enttype %in% enttype_hdl)
  df_output_hdl[[p]] <- df
  
  # LDL cholesterol
  df <- filter(alpha, enttype %in% enttype_ldl)
  df_output_ldl[[p]] <- df
  
  # Triglycerides
  df <- filter(alpha, enttype %in% enttype_trig)
  df_output_trig[[p]] <- df
  
  # Hba1c
  df <- filter(alpha, enttype %in% enttype_hba1c)
  df_output_hba1c[[p]] <- df
  
  # Haemoglobin
  df <- filter(alpha, enttype %in% enttype_haem)
  df_output_haem[[p]] <- df
  
  # Calcium
  df <- filter(alpha, enttype %in% enttype_calc)
  df_output_calc[[p]] <- df
  
  # ALT
  df <- filter(alpha, enttype %in% enttype_alt)
  df_output_alt[[p]] <- df
  
  # ALP
  df <- filter(alpha, enttype %in% enttype_alp)
  df_output_alp[[p]] <- df
  
  # Bilirubin
  df <- filter(alpha, enttype %in% enttype_bil)
  df_output_bil[[p]] <- df
  
  # Phosphate
  df <- filter(alpha, enttype %in% enttype_phos)
  df_output_phos[[p]] <- df
  
  # Serum Albumin
  df <- filter(alpha, enttype %in% enttype_alb)
  df_output_alb[[p]] <- df
  
}

# combine & save

setwd(outputDir)

df_chol <- do.call(rbind, df_output_chol)
print(nrow(df_chol))
save(df_chol, file = "df_chol.Rdata")

df_chol <- do.call(rbind, df_output_chol)
print(nrow(df_chol))
save(df_chol, file = "df_serum_chol.Rdata")

df_hdl <- do.call(rbind, df_output_hdl)
print(nrow(df_hdl))
save(df_hdl, file = "df_hdl.Rdata")

df_ldl <- do.call(rbind, df_output_ldl)
print(nrow(df_ldl))
save(df_ldl, file = "df_ldl.Rdata")

df_trig <- do.call(rbind, df_output_trig)
print(nrow(df_trig))
save(df_trig, file = "df_trig.Rdata")

df_hba1c <- do.call(rbind, df_output_hba1c)
print(nrow(df_hba1c))
save(df_hba1c, file = "df_hba1c.Rdata")

df_haem <- do.call(rbind, df_output_haem)
print(nrow(df_haem))
save(df_haem, file = "df_haem.Rdata")

df_calc <- do.call(rbind, df_output_calc)
print(nrow(df_calc))
save(df_calc, file = "df_calc.Rdata")

df_alt <- do.call(rbind, df_output_alt)
print(nrow(df_alt))
save(df_alt, file = "df_alt.Rdata")

df_alp <- do.call(rbind, df_output_alp)
print(nrow(df_alp))
save(df_alp, file = "df_alp.Rdata")

df_bil <- do.call(rbind, df_output_bil)
print(nrow(df_bil))
save(df_bil, file = "df_bil.Rdata")

df_phos <- do.call(rbind, df_output_phos)
print(nrow(df_phos))
save(df_phos, file = "df_phos.Rdata")

df_alb <- do.call(rbind, df_output_alb)
print(nrow(df_alb))
save(df_alb, file = "df_serum_alb.Rdata")


######################################################################
### QRISK variables: Clinical + additional
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)

sourceDir_codelist <- "L:\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "D:\\CKD\\source\\"
sourceDir_R <- "L:\\Monitoring CKD CHF\\Data\\derived CPRD data\\"
outputDir <- "L:\\Monitoring CKD CHF\\Data\\derived CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### entity

setwd(sourceDir_R)
load("entity.Rdata")

### entity types

source("L:\\Monitoring CKD CHF\\Data\\code\\functions.R")
colname <- "medcode"

# Blood pressure
enttype_bp <- subset(entity, description == "Blood pressure")$enttype # 1
df_output_bp <- NULL

# BMI
enttype_height <- filter(entity, description == "Height")$enttype # entity 14
enttype_weight <- filter(entity, description == "Weight")$enttype # entity 13; includes BMI
df_output_height <- NULL
df_output_weight <- NULL

# Chol/HDL ratio
enttype_chol <- subset(entity, description == "HDL/LDL ratio")$enttype # 338
df_output_chol <- NULL

# smoking
medcodes_smok <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "smoking_final.csv", colname = colname)
df_output_smok <- NULL

### extract data

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  ### Clinical file ###
  
  setwd(paste(sourceDir_CPRD, "Clinical", sep = "\\"))
  f <- paste("p", p, "_Clinical_cut.dta", sep = "")
  df_clinical <- read.dta13(f)
  head(df_clinical)
  df_clinical <- select(df_clinical, patid, eventdate, adid, medcode, enttype)
  df_clinical <- mutate(df_clinical, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  ### Additional file ###
  
  setwd(paste(sourceDir_CPRD, "Additional", sep = "\\"))
  f <- paste("p", p, "_Additional_cut.dta", sep = "")
  df_add <- read.dta13(f)
  
  ### merge and clean ###
  
  alpha <- merge(df_clinical, df_add)
  
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove duplicates
  alpha <- distinct(alpha)
  
  # SBP
  df <- filter(alpha, enttype %in% enttype_bp)
  df_output_bp[[p]] <- df
  print(paste("BP:", nrow(df)))
  
  # Height
  df <- filter(alpha, enttype %in% enttype_height)
  df_output_height[[p]] <- df
  print(paste("Height:", nrow(df)))
  
  # Weight
  df <- filter(alpha, enttype %in% enttype_weight)
  df_output_weight[[p]] <- df
  print(paste("Weight:", nrow(df)))
  
  # smoking
  df <- filter(alpha, medcode %in%  medcodes_smok)
  df_output_smok[[p]] <- df
  print(paste("smoking:", nrow(df)))
}

# combine
df_bp <- do.call(rbind, df_output_bp)
print(nrow(df_bp))
df_height <- do.call(rbind, df_output_height)
print(nrow(df_height))
df_weight <- do.call(rbind, df_output_weight)
print(nrow(df_weight))
df_smok <- do.call(rbind, df_output_smok)
print(nrow(df_smok))

# save
setwd(outputDir)
save(df_bp, file = "df_bp.Rdata")
save(df_height, file = "df_height.Rdata")
save(df_weight, file = "df_weight.Rdata")
save(df_smok, file = "df_smok.Rdata")

######################################################################
### QRISK variables: prescriptions
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\derived CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### medcodes ###

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")
colname <- "prodcode"

# statins
prodcodes_statins <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "statins_final.csv", colname = "prodcode")
df_output_statins <- NULL

# insulins
prodcodes_ins <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "insulins_final.csv", colname = "prodcode")
df_output_ins <- NULL

# anti-hypertensives
prodcodes_antihyp <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "antihypertensives_final.csv", colname = "prodcode")
df_output_antihyp <- NULL

# steroids
prodcodes_steroids <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "steroids_final.csv", colname = "prodcode")
df_output_steroids <- NULL

# antipsychotic
prodcodes_antipsych <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "antipsych_final.csv", colname = "prodcode")
df_output_antipsych <- NULL

# erectile dysfunction
prodcodes_erectile_Tx <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "erectile_Tx_final.csv", colname = "prodcode")
df_output_erectile_Tx <- NULL

# ACEI
prodcodes_acei <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "acei_final.csv", colname = "prodcode")
df_output_acei <- NULL

# ARB
prodcodes_arb <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "arb_final.csv", colname = "prodcode")
df_output_arb <- NULL

# diuretics
prodcodes_diur <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "diuretics_final.csv", colname = "prodcode")
df_output_diur <- NULL

# NSAIDS
prodcodes_nsaid <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "nsaid_final.csv", colname = "prodcode")
df_output_nsaid <- NULL

# Phosphate-binder medications
prodcodes_pb <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "pb_final.csv", colname = "prodcode")
df_output_pb <- NULL

# Anti-platelets and anti-coagulant medications
prodcodes_apc <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "Antiplatelet_Anticoagulants.csv", colname = "prodcode")
df_output_apc <- NULL

### extract data

# local function for keeping only patients with at least two prescriptions
.extract_df <- function(alpha, prodcodes_t) {
  df <- filter(alpha, prodcode %in% prodcodes_t)
  df <- group_by(df, patid)
  df_2 <- summarise(df, n = length(eventdate))
  df_2 <- filter(df_2, n >= 2)
  id <- df_2$patid
  return(filter(df, patid %in% id))
}

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  # load the practice dataset: Therapy
  setwd(paste(sourceDir_CPRD, "Therapy", sep = "\\"))
  f <- paste("p", p, "_Therapy_cut.dta", sep = "")
  df_therapy <- read.dta13(f)
  df_therapy <- select(df_therapy, patid, prodcode, eventdate)
  df_therapy <- mutate(df_therapy, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  ### initialise output and do basic cleaning
  
  alpha <- df_therapy
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove duplicates
  alpha <- distinct(alpha)
  
  # statins
  #df <- .extract_df(alpha, prodcodes_statins)
  #df_output_statins[[p]] <- df
  #print(paste("statins:", nrow(df)))
  
  # insulins
  #df <- .extract_df(alpha, prodcodes_ins)
  #df_output_ins[[p]] <- df
  #print(paste("insulins:", nrow(df)))
  
  # anti-hypertensives
  #df <- .extract_df(alpha, prodcodes_antihyp)
  #df_output_antihyp[[p]] <- df
  #print(paste("antihyp:", nrow(df)))
  
  # steroids
  #df <- .extract_df(alpha, prodcodes_steroids)
  #df_output_steroids[[p]] <- df
  #print(paste("steroids:", nrow(df)))
  
  # antipsych
  #df <- .extract_df(alpha, prodcodes_antipsych)
  #df_output_antipsych[[p]] <- df
  #print(paste("antipsych:", nrow(df)))
  
  # erectile Tx
  #df <- .extract_df(alpha, prodcodes_erectile_Tx)
  #df_output_erectile_Tx[[p]] <- df
  #print(paste("erectile Tx:", nrow(df)))
  
  # ACEI
  #df <- .extract_df(alpha, prodcodes_acei)
  #df_output_acei[[p]] <- df
  
  # ARB
  #df <- .extract_df(alpha, prodcodes_arb)
  #df_output_arb[[p]] <- df
  
  # diuretics
  #df <- .extract_df(alpha, prodcodes_diur)
  #df_output_diur[[p]] <- df
  
  # NSAIDs
  #df <- .extract_df(alpha, prodcodes_nsaid)
  #df_output_nsaid[[p]] <- df
  
  # Phosphate binders
  #df <- .extract_df(alpha, prodcodes_pb)
  #df_output_pb[[p]] <- df
  
  # Antiplatelets or anti-coagulants
  df <- .extract_df(alpha, prodcodes_apc)
  df_output_apc[[p]] <- df
  print(paste("apc:", nrow(df)))
  
}

# combine & save

setwd(outputDir)

df_statins <- do.call(rbind, df_output_statins)
print(nrow(df_statins))
#save(df_statins, file = "df_statins.Rdata")

df_ins <- do.call(rbind, df_output_ins)
print(nrow(df_ins))
#save(df_ins, file = "df_ins.Rdata")

df_antihyp <- do.call(rbind, df_output_antihyp)
print(nrow(df_antihyp))
#save(df_antihyp, file = "df_antihyp.Rdata")

df_steroids <- do.call(rbind, df_output_steroids)
print(nrow(df_steroids))
#save(df_steroids, file = "df_steroids.Rdata")

df_antipsych <- do.call(rbind, df_output_antipsych)
print(nrow(df_antipsych))
#save(df_antipsych, file = "df_antipsych.Rdata")

df_erectile_Tx <- do.call(rbind, df_output_erectile_Tx)
print(nrow(df_erectile_Tx))
#save(df_erectile_Tx, file = "df_erectile_Tx.Rdata")

df_acei <- do.call(rbind, df_output_acei)
print(nrow(df_acei))
#save(df_acei, file = "df_acei.Rdata")

df_arb <- do.call(rbind, df_output_arb)
print(nrow(df_arb))
#save(df_arb, file = "df_arb.Rdata")

df_diur <- do.call(rbind, df_output_diur)
print(nrow(df_diur))
#save(df_diur, file = "df_diur.Rdata")

df_nsaid <- do.call(rbind, df_output_nsaid)
print(nrow(df_nsaid))
#save(df_nsaid, file = "df_nsaid.Rdata")

df_pb <- do.call(rbind, df_output_pb)
print(nrow(df_pb))
#save(df_pb, file = "df_pb.Rdata")

df_apc <- do.call(rbind, df_output_apc)
print(nrow(df_apc))
save(df_apc, file = "df_apc_CLS.Rdata")

######################################################################
### Additional extraction - to be combined with the above text
######################################################################


######################################################################
### Clinical
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)
library(data.table)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### medcodes ###

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")
colname <- "medcode"

## transplant
#medcodes_transplant <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "renal_transplant_final.csv", colname = colname)
#df_output_transplant <- NULL

# pregnancy
medcodes_pregnancy <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "pregnancy_final.csv", colname = colname)
df_output_pregnancy <- NULL

### extract data

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  # load the practice dataset: clinical
  setwd(paste(sourceDir_CPRD, "Clinical", sep = "\\"))
  f <- paste("p", p, "_Clinical_cut.dta", sep = "")
  df_clinical <- read.dta13(f)
  df_clinical <- select(df_clinical, patid, eventdate, medcode)
  df_clinical <- mutate(df_clinical, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  # load the practice dataset: referral
  setwd(paste(sourceDir_CPRD, "Referral", sep = "\\"))
  f <- paste("p", p, "_Referral_cut.dta", sep = "")
  df_ref <- read.dta13(f)
  df_ref <- select(df_ref, patid, eventdate, medcode)
  df_ref <- mutate(df_ref, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  ### initialise output and do basic cleaning
  
  alpha <- rbind(df_clinical, df_ref)
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove duplicates
  alpha <- distinct(alpha)
  
  ## transplant
  #df <- filter(alpha, medcode %in% medcodes_transplant)
  #df_output_transplant[[p]] <- df
  #print(paste("Transplant:", nrow(df)))
  
  # pregnancy
  df <- filter(alpha, medcode %in%  medcodes_pregnancy)
  df_output_pregnancy[[p]] <- df
  print(paste("Pregnancy:", nrow(df)))
  
}

# combine and save

setwd(outputDir)

#df_transplant <- do.call(rbind, df_output_transplant)
#print(nrow(df_transplant))
#save(df_transplant, file = "df_transplant.Rdata")

df_pregnancy <- do.call(rbind, df_output_pregnancy)
print(nrow(df_pregnancy))
save(df_pregnancy, file = "df_pregnancy.Rdata")

######################################################################
### Clinical  + additional
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)
library(data.table)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### medcodes ###

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")

# FH CHD
df_output_fh_chd <- NULL
enttype_fh_chd <- 87

### extract data

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  ### Clinical file ###
  
  setwd(paste(sourceDir_CPRD, "Clinical", sep = "\\"))
  f <- paste("p", p, "_Clinical_cut.dta", sep = "")
  df_clinical <- data.table(read.dta13(f))
  df_clinical
  df_clinical <- df_clinical[, .(patid, eventdate, medcode, enttype, adid)]
  df_clinical[, eventdate := as.Date(eventdate, format = "%d/%m/%Y")]
  df_clinical <- df_clinical[!is.na(eventdate)]
  df_clinical <- df_clinical[enttype == enttype_fh_chd]
  
  ### Additional file ###
  
  setwd(paste(sourceDir_CPRD, "Additional", sep = "\\"))
  f <- paste("p", p, "_Additional_cut.dta", sep = "")
  df_add <- data.table(read.dta13(f))
  
  ### merge and clean ###
  
  alpha <- merge(df_clinical, df_add, by = c("patid", "enttype", "adid"))
  alpha <- distinct(alpha)
  
  # pregnancy
  df <- alpha
  df_output_fh_chd[[p]] <- df
  print(paste("FH CHD:", nrow(df)))
  
}

# combine and save

setwd(outputDir)

#df_transplant <- do.call(rbind, df_output_transplant)
#print(nrow(df_transplant))
#save(df_transplant, file = "df_transplant.Rdata")

df_fh_chd <- do.call(rbind, df_output_fh_chd)
print(nrow(df_fh_chd))
save(df_fh_chd, file = "df_chd_2.Rdata")

######################################################################
### Albuminuria codes
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)
library(data.table)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### medcodes ###

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")

# albuminuria codes
setwd(sourceDir_codelist)
codes_t <- read.csv("albuminuria_all_BF_final.csv")
medcodes_alb <- codes_t$medcode[codes_t$included == 1]
df_output_alb_clin <- NULL
df_output_alb_test <- NULL

### extract data

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  ### Clinical file ###
  
  setwd(paste(sourceDir_CPRD, "Clinical", sep = "\\"))
  f <- paste("p", p, "_Clinical_cut.dta", sep = "")
  df_clinical <- data.table(read.dta13(f))
  df_clinical
  df_clinical <- df_clinical[, .(patid, eventdate, medcode, enttype, adid)]
  df_clinical[, eventdate := as.Date(eventdate, format = "%d/%m/%Y")]
  df_clinical <- df_clinical[!is.na(eventdate)]
  df_clinical <- df_clinical[medcode %in% medcodes_alb]
  
  ### Additional file ###
  
  setwd(paste(sourceDir_CPRD, "Additional", sep = "\\"))
  f <- paste("p", p, "_Additional_cut.dta", sep = "")
  df_add <- data.table(read.dta13(f))
  
  ### merge and clean ###
  
  alpha <- merge(df_clinical, df_add, by = c("patid", "enttype", "adid"), all.x = TRUE)
  alpha <- distinct(alpha)
  
  # save
  df_output_alb_clin[[p]] <- alpha
  print(paste("Alb (diagnosis):", nrow(alpha)))
  
  ### Test file ###

  setwd(paste(sourceDir_CPRD, "Test", sep = "\\"))
  
  # load the practice dataset: test
  f <- paste("p", p, "_Test_cut.dta", sep = "")
  df_test <- data.table(read.dta13(f))
  df_test <- df_test[medcode %in% medcodes_alb]
  df_test <- df_test[, .(patid, eventdate, medcode, enttype, data1, data2, data3, data4)]
  df_test[, eventdate := as.Date(eventdate, format = "%d/%m/%Y")]
  # 
  # initialise output and do basic cleaning
  # remove tests with no date
  alpha <- df_test[!is.na(eventdate)]
  ## remove tests with empty value <- NOT HERE, as tests may have different formats!
  #alpha <- filter(alpha, !is.na(data2))
  # remove duplicates
  alpha <- distinct(alpha)
  
  # save
  df_output_alb_test[[p]] <- alpha
  print(paste("Alb (test):", nrow(alpha)))
  
}

# combine and save

setwd(outputDir)

df_alb_clin <- do.call(rbind, df_output_alb_clin)
print(nrow(df_alb_clin))
save(df_alb_clin, file = "df_alb_clin.Rdata")

df_alb_test <- do.call(rbind, df_output_alb_test)
print(nrow(df_alb_test))
save(df_alb_test, file = "df_alb_test.Rdata")

######################################################################
### ACR entity: test file
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(dplyr)
library(data.table)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### entity

setwd(sourceDir_R)
load("entity.Rdata")

### entity types

# ACR
enttype_acr <- subset(entity, description == "Albumin creatinine ratio")$enttype # 469
df_output_acr <- NULL

### extract data
i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  ### Test file ###
  
  setwd(paste(sourceDir_CPRD, "Test", sep = "\\"))
  
  # load the practice dataset: test
  f <- paste("p", p, "_Test_cut.dta", sep = "")
  df_test <- data.table(read.dta13(f))
  df_test <- df_test[, .(patid, eventdate, medcode, enttype, data1, data2, data3)]
  df_test[, eventdate := as.Date(eventdate, format = "%d/%m/%Y")]
  
  # initialise output and do basic cleaning
  alpha <- df_test
  # remove tests with no date
  alpha <- alpha[!is.na(eventdate)]
  ## remove tests with empty value <- NOT HERE, as tests may have different formats!
  #alpha <- filter(alpha, !is.na(data2))
  # remove duplicates
  alpha <- distinct(alpha)
  
  # ACR
  df <- alpha[enttype %in% enttype_acr]
  df_output_acr[[p]] <- df
  print(paste("ACR:", nrow(df)))
  
}

# combine & save

setwd(outputDir)

df_acr <- do.call(rbind, df_output_acr)
print(nrow(df_acr))
save(df_acr, file = "df_acr_entity.Rdata")

######################################################################
### QRISK variables: prescriptions
######################################################################

### load the data ###

rm(list = ls())

library(readstata13)
library(plyr)
library(dplyr)
library(gdata)
#library(data.table)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_CPRD <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\CKD\\source\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"

### list of eligible practices

setwd(sourceDir_R)
load("link_pracid.Rdata")

### medcodes ###

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")
colname <- "prodcode"

# beta-blockers
prodcodes_beta <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "beta_final.csv", colname = "prodcode")
df_output_beta <- NULL

# ccb
prodcodes_ccb <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = "ccb_final.csv", colname = "prodcode")
df_output_ccb <- NULL

### extract data

# local function for keeping only patients with at least two prescriptions
.extract_df <- function(alpha, prodcodes_t) {
  df <- filter(alpha, prodcode %in% prodcodes_t)
  df <- group_by(df, patid)
  df_2 <- summarise(df, n = length(eventdate))
  df_2 <- filter(df_2, n >= 2)
  id <- df_2$patid
  return(filter(df, patid %in% id))
}

p <- link_pracid[1]

i <- 0
for (p in link_pracid){
  i <- i + 1
  print(paste("Practice = ", p, sep = ""))
  print(i)
  
  # load the practice dataset: Therapy
  setwd(paste(sourceDir_CPRD, "Therapy", sep = "\\"))
  f <- paste("p", p, "_Therapy_cut.dta", sep = "")
  df_therapy <- read.dta13(f)
  df_therapy <- select(df_therapy, patid, prodcode, eventdate)
  df_therapy <- mutate(df_therapy, eventdate = as.Date(eventdate, format = "%d/%m/%Y"))
  
  ### initialise output and do basic cleaning
  
  alpha <- df_therapy
  # remove tests with no date
  alpha <- filter(alpha, !is.na(eventdate))
  # remove duplicates
  alpha <- distinct(alpha)
  
  # beta-blockers
  df <- .extract_df(alpha, prodcodes_beta)
  df_output_beta[[p]] <- df
  print(paste("beta-blockers:", nrow(df)))
  
  # CCB
  df <- .extract_df(alpha, prodcodes_ccb)
  df_output_ccb[[p]] <- df
  print(paste("CCB:", nrow(df)))
}

# combine & save

setwd(outputDir)

df_beta <- do.call(rbind, df_output_beta)
print(nrow(df_beta))
save(df_beta, file = "df_beta.Rdata")

df_ccb <- do.call(rbind, df_output_ccb)
print(nrow(df_ccb))
save(df_ccb, file = "df_ccb.Rdata")