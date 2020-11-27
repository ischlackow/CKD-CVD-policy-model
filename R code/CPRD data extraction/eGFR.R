######################################################################
### Setup
######################################################################

rm(list = ls())

library(readstata13)
library(dplyr)
library(gdata)
library(data.table)

sourceDir_lookup <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\TXTFILES\\"
sourceDir_R_raw <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"

######################################################################
### Load important files
######################################################################

### patient file

setwd(sourceDir_R)
df_pat <- data.table(get(load("df_pat.Rdata")))
df_pat <- df_pat[, .(patid, sex, dob)]

### CPRD lookups

setwd(sourceDir_lookup)
OPR <- data.table(read.delim("OPR.txt", stringsAsFactors = FALSE))
SUM <- data.table(read.delim("SUM.txt", stringsAsFactors = FALSE))

# creatinine tests
setwd(sourceDir_R_raw)
df_scr <- data.table(get(load("df_scr.Rdata")))

######################################################################
### Clean creatinine tests
######################################################################

### DO NOT filter by censoring dates

### operator
df_scr <- rename(df_scr, Code = data1)
df_scr <- merge(df_scr, OPR, by = "Code", all.x = TRUE)
df_scr[, Code := NULL]
df_scr[, Operator := drop.levels(Operator)]
table(df_scr$Operator, exclude = NULL)
# leave only = for the moment; this is a majority anyway
df_scr <- df_scr[Operator == "="]
df_scr[, Operator := NULL]

### units
df_scr <- rename(df_scr, Code = data3)
df_scr <- merge(df_scr, SUM, by = "Code")
df_scr[, Code := NULL]
df_scr[, Specimen.Unit.Of.Measure := drop.levels(df_scr$Specimen.Unit.Of.Measure)]
table(df_scr$Specimen.Unit.Of.Measure, exclude = NULL)
# leave only those that indicate mg/L or mol/L
df_scr <- df_scr[Specimen.Unit.Of.Measure %in% c("g/L", "mg/dL", "mg/L", "mmol/L", "mol/L", "ng/mL", "nmol/L", "pmol/L", "ug/L", "umol/L")]

# check plausibility and convert everything into mg/dL

df_scr[, scr.mg.dl := NaN]

# g/L
df_scr[Specimen.Unit.Of.Measure == "g/L", scr.mg.dl :=  data2 * 100]
quantile(df_scr[Specimen.Unit.Of.Measure == "g/L", scr.mg.dl]) # consistently too large; remove
df_scr <- df_scr[Specimen.Unit.Of.Measure != "g/L"]

# mg/dL - leave as is
df_scr[Specimen.Unit.Of.Measure == "mg/dL", scr.mg.dl := data2]
quantile(df_scr[Specimen.Unit.Of.Measure == "mg/dL", scr.mg.dl])

# mg/L
df_scr[Specimen.Unit.Of.Measure == "mg/L", scr.mg.dl := data2 * 0.1]
quantile(df_scr[Specimen.Unit.Of.Measure == "mg/L", scr.mg.dl])

# mmol/L
df_scr[Specimen.Unit.Of.Measure == "mmol/L", scr.mg.dl := data2 / 0.08842]
quantile(df_scr[Specimen.Unit.Of.Measure == "mmol/L", scr.mg.dl])
quantile(df_scr[Specimen.Unit.Of.Measure == "mmol/L", scr.mg.dl], probs = c(0.01, 0.05, 0.1)) # even the bottom 1% is high; remove
df_scr <- df_scr[Specimen.Unit.Of.Measure != "mmol/L"]

# mol/L
df_scr[Specimen.Unit.Of.Measure == "mol/L", scr.mg.dl := data2 * 1000 / 0.08842]
quantile(df_scr[Specimen.Unit.Of.Measure == "mol/L", scr.mg.dl]) # consistently too large; remove
df_scr <- df_scr[Specimen.Unit.Of.Measure != "mol/L"]

# ng/mL
df_scr[Specimen.Unit.Of.Measure == "ng/mL", scr.mg.dl := data2 * 0.0001]
quantile(df_scr[Specimen.Unit.Of.Measure == "ng/mL", scr.mg.dl]) # seems on the lower side; will probably be removed at the later filtering, so leave for now

# nmol/L
df_scr[Specimen.Unit.Of.Measure == "nmol/L", scr.mg.dl := data2 * 0.000001 / 88.4]
quantile(df_scr[Specimen.Unit.Of.Measure == "nmol/L", scr.mg.dl]) # too low; remove
df_scr <- df_scr[Specimen.Unit.Of.Measure != "nmol/L"]

# pmol/L
df_scr[Specimen.Unit.Of.Measure == "pmol/L", scr.mg.dl := data2 * 0.000000001 / 0.08842]
quantile(df_scr[Specimen.Unit.Of.Measure == "pmol/L", scr.mg.dl]) # too low; remove
df_scr <- df_scr[Specimen.Unit.Of.Measure != "pmol/L"]

# ug/L
df_scr[Specimen.Unit.Of.Measure == "ug/L", scr.mg.dl := data2 * 0.001 * 0.1]
quantile(df_scr[Specimen.Unit.Of.Measure == "ug/L", scr.mg.dl]) # possibly too low; leave until the next filtering

# umol/L
df_scr[Specimen.Unit.Of.Measure == "umol/L", scr.mg.dl := data2 / 88.4]
quantile(df_scr[Specimen.Unit.Of.Measure == "umol/L", scr.mg.dl]) # possibly too low; leave until the next filtering

# check that all values have been converted and remove those <10 umol/L (as per DL's advice)
range(df_scr$scr.mg.dl)
df_scr <- df_scr[scr.mg.dl >= 10/88.4]
df_scr[, Specimen.Unit.Of.Measure := NULL]
df_scr[, data2 := NULL]
quantile(df_scr$scr.mg.dl, probs = c(0.8, 0.9, 0.99)) # so only a very few abnormally high values
range(df_scr$scr.mg.dl)

### average tests taken on the same day ###

df_scr <- df_scr[, lapply(.SD, mean), by = .(patid, eventdate), .SDcols = c("scr.mg.dl")]
range(df_scr$scr.mg.dl)

######################################################################
### Calculate eGFR
######################################################################

### add information on whether the patient is black ###

setwd(sourceDir_R)
ids_black <- get(load("ids_black.Rdata"))
df_scr[, black := as.numeric(patid %in% ids_black)]
table(df_scr$black, exclude = NULL)

### add sex and age ###

df_scr <- merge(df_scr, df_pat, by = "patid")
df_scr[, age := as.numeric(difftime(eventdate, dob, units = "days") / 365.25)]
range(df_scr$age)
# leave only tests taken at or after the age of 18
df_scr <- df_scr[age >= 18]
df_scr

# ckd-epi: p14 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2763564/pdf/nihms132246.pdf
df_scr[, k := ifelse(sex == "F", 0.7, 0.9)]
df_scr[, a := ifelse(sex == "F", -0.329, -0.411)]
df_scr[, egfr := 141 * pmin(scr.mg.dl/k, 1)^a * pmax(scr.mg.dl/k, 1)^(-1.209) * (0.993^(age))]
df_scr[sex == "F", egfr := 1.018 * egfr]
df_scr[black == 1, egfr := 1.159 * egfr]

# stage
df_scr[, ckd := cut(egfr, breaks = c(0, 15, 30, 45, 60, 90, Inf), labels = c("ckd5", "ckd4", "ckd3b", "ckd3a", "ckd2", "ckd1"))]
df_scr[, ckd_num := as.numeric(factor(df_scr$ckd, levels(df_scr$ckd)[6 : 1]))]
table(df_scr$ckd, df_scr$ckd_num)

# 6 = ckd5
# 5 = ckd4
# 4 = ckd3b
# 3 = ckd3a
# 2 = ckd2
# 1 = ckd1 or no ckd

# check
alpha <- df_scr[, lapply(.SD, mean), by = ckd_num, .SDcols = c("egfr")]
alpha

### leave only patients with at least two egfr measurements of CKD stage 2-6

alpha <- df_scr
length(unique(alpha$patid))
alpha[, has_ckd := as.numeric(ckd_num >= 2)]
table(alpha$has_ckd)
alpha <- alpha[, lapply(.SD, sum), by = .(patid), .SDcols = c("has_ckd")]
alpha[has_ckd>=1]
alpha <- alpha[has_ckd >= 2]
id <- alpha$patid
df_scr <- df_scr[patid %in% id]

### save ###

df_egfr <- df_scr[, .(patid, eventdate, egfr, ckd, ckd_num)]
setwd(outputDir)
save(df_egfr, file = "df_egfr.Rdata")
