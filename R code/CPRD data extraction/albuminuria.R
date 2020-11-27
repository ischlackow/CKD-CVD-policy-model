# Note: definition of albuminuria/proteinuria stages is based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4089693/table/tbl7/

######################################################################
### Setup
######################################################################

rm(list = ls())

library(readstata13)
library(dplyr)
library(gdata)
library(data.table)

sourceDir_lookup <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\TXTFILES\\"
sourceDir_codelists <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_R_raw <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\albuminuria_functions.R")

######################################################################
### Files
######################################################################

# CPRD lookups
setwd(sourceDir_lookup)
OPR <- read.delim("OPR.txt", stringsAsFactors = FALSE)
SUM <- read.delim("SUM.txt", stringsAsFactors = FALSE)
TQU <- read.delim("TQU.txt", stringsAsFactors = FALSE)

# entity
setwd(sourceDir_R)
entity <- data.table(get(load("entity.Rdata")))

######################################################################
### Output dataset
######################################################################

output <- NULL

######################################################################
### Numeric values: Albumin-to-creatinine ratio
######################################################################

### Scenario 1: stand-alone ACR test ###

# initial cleaning
df <- .clean_df(filename = "df_acr.Rdata")

# units & conversion

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# check when data2 = 0
test <- df[data2 == 0]
table(test$Specimen.Unit.Of.Measure, exclude = NULL) # keep those that look vaguely correct

# remove only units vaguely similar to mg/mmol or mg/g
# KEEP units that have (creat) in them
df <- df[(Specimen.Unit.Of.Measure %in% c("g/mol", "mg/mmol", "mg/mmol(creat)", "ug/mmol(creat)")) | 
           (data2 == 0 & Specimen.Unit.Of.Measure %in% c("g/mol", "mg/mmol", "mg/mmol(creat)", "No Data Entered", "ratio"))]
#nrow(df[data2 == 0]) - nrow(test)
# convert all to mg/mmol
# in fact, g/mol is the same as mg/mmol, and the same true for data2 = 0
# only need to convert ug/mmol
df[, acr.mg.mmol := data2]
df[Specimen.Unit.Of.Measure == "ug/mmol(creat)", acr.mg.mmol := data2 * 0.001]

quantile(df[Specimen.Unit.Of.Measure %in% c("g/mol"), acr.mg.mmol])
quantile(df[Specimen.Unit.Of.Measure %in% c("mg/mmol"), acr.mg.mmol]) # possibly a few abnormal values, but nothing consistently wrong
quantile(df[Specimen.Unit.Of.Measure %in% c("mg/mmol(creat)"), acr.mg.mmol])
quantile(df[Specimen.Unit.Of.Measure %in% c("ug/mmol(creat)"), acr.mg.mmol]) # probably too low but leave as is

# interpretable values

table(df$Operator, exclude = NULL)
lb <- 3
ub <- 30
df <- .filter_lb_up(df, colname = "acr.mg.mmol", lb = lb, ub = ub)

# plausible values

range(df$acr.mg.mmol)
subset(df, acr.mg.mmol >= 10000) 
# remove >50000
df <- df[acr.mg.mmol <= 50000]

# final clean

df <- .clean_df_post(df = df, colname = "acr.mg.mmol", mean_per_day = FALSE)
df[, Operator := NULL]
df <- rename(df, val = acr.mg.mmol)

# assign stage 

df <- .assign_stage_cts(df, lb, ub)
table(df$stage, exclude = NULL)

# add to the output

df[, test := "acr.mg.mmol"]
df[, lab := "acr, standalone"]
output <- rbind(output, df)

### Scenario 2: separate urine albumin & urine creatinine tests ###

### urine albumin concentrations

df <- .clean_df(filename = "df_alb.Rdata")

# units & conversion

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# check when data2 = 0
test <- df[data2 == 0]
table(test$Specimen.Unit.Of.Measure, exclude = NULL) # all apart from %

# leave only units vaguely look like mg/L
df <- df[(Specimen.Unit.Of.Measure %in% c("g/dL", "g/L", "mg/dL", "mg/L",  "mg/m3", "ng/L")) |
           (data2 == 0 & Specimen.Unit.Of.Measure %in% c("g/L", "mg/L", "No Data Entered"))]
#nrow(df[data2 == 0]) - nrow(test)
df[, alb.mg.L := NaN]

df[Specimen.Unit.Of.Measure %in% c("g/dL"), alb.mg.L := data2 * 10 * 1000]
quantile(df$alb.mg.L[df$Specimen.Unit.Of.Measure %in% c("g/dL")]) # consistently too high; remove
df <- df[Specimen.Unit.Of.Measure != "g/dL"]

df[Specimen.Unit.Of.Measure %in% c("g/L"), alb.mg.L := data2 * 1000]
quantile(df$alb.mg.L[df$Specimen.Unit.Of.Measure %in% c("g/L")]) # some values too high, but some could be plausible, so leave

df[Specimen.Unit.Of.Measure %in% c("mg/dL"), alb.mg.L := data2 * 10]
quantile(df$alb.mg.L[df$Specimen.Unit.Of.Measure %in% c("mg/dL")])

df[Specimen.Unit.Of.Measure %in% c("mg/L"), alb.mg.L := data2]
quantile(df$alb.mg.L[df$Specimen.Unit.Of.Measure %in% c("mg/L")])

df[Specimen.Unit.Of.Measure %in% c("mg/m3"), alb.mg.L := data2 * 0.001]
quantile(df$alb.mg.L[df$Specimen.Unit.Of.Measure %in% c("mg/m3")]) # consistently low; remove
df <- df[Specimen.Unit.Of.Measure != "mg/m3"]

df[Specimen.Unit.Of.Measure %in% c("ng/L"), alb.mg.L := data2 * 0.000001]
quantile(df$alb.mg.L[df$Specimen.Unit.Of.Measure %in% c("ng/L")]) # consistently low; remove
df <- df[Specimen.Unit.Of.Measure != "ng/L"]

df[data2 == 0, alb.mg.L := 0]

# plausible values

range(df$alb.mg.L)
subset(df, alb.mg.L >= 50000) # over 60,000 seems to be consistently too high
subset(df, patid %in% subset(df, alb.mg.L >= 60000)$patid) # seems a misprint; remove
df <- df[alb.mg.L <= 60000]

# interpretable values

table(df$Operator, exclude = NULL)
# leave all as is; filter at the ratio stage

# final clean & save

df_alb <- .clean_df_post(df, colname = "alb.mg.L", mean = FALSE)

### urine creatinine concentrations

## NB: it seems that urine creatinine concentrations are meant to be measured in mmol/day. Hoowever, the vast majority of the tests are in  mmol/L & umol/L, 
# it is not clear whether these are maybe serum creatinine tests?
# for the moment leave mmol/L & umol/L; check how many additional tests this gives.
# See also http://www.mps.com.au/media/1648771/urine%20insight%20final%20low%20res.pdf for the cutoff value of 1.8 mmol/L, albeit in a different context, but still a urine creatinine test

df <- .clean_df("df_cr.Rdata")

# remove zeroes, as they are going into the denominator

df <- df[data2 > 0]

# units & conversion

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# leave only units vaguely look like mmol/L or g/L
# cannot leave No Data Entered as only non-zero values are considered
df <- df[Specimen.Unit.Of.Measure %in% c("g/L", "mg/dL", "mg/L", "mg/m3", "mmol/L", "mol/L", "nmol/L", "umol/L")]

df[, cr.mmol.L := NaN]

df[Specimen.Unit.Of.Measure %in% c("g/L"), cr.mmol.L := data2 * 100 * 88.4]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("g/L")])

df[Specimen.Unit.Of.Measure %in% c("mg/dL"), cr.mmol.L := data2 * 88.4]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("mg/dL")])

df[Specimen.Unit.Of.Measure %in% c("mg/L"), cr.mmol.L := data2 * 0.1 * 88.4]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("mg/L")])

df[Specimen.Unit.Of.Measure %in% c("mg/m3"), cr.mmol.L := data2 * 0.001 * 0.1 * 88.4]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("mg/m3")])

df[Specimen.Unit.Of.Measure %in% c("mmol/L"), cr.mmol.L := data2]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("mmol/L")]) # some very high values
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("mmol/L")], probs = c(0.8, 0.9, 0.99))

df[Specimen.Unit.Of.Measure %in% c("mol/L"), cr.mmol.L := data2 * 1000]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("mol/L")]) # some very high values

df[Specimen.Unit.Of.Measure %in% c("nmol/L"), cr.mmol.L := data2 * 0.000001]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("nmol/L")]) # all too small; remove
df <- df[Specimen.Unit.Of.Measure != "nmol/L"]

df[Specimen.Unit.Of.Measure %in% c("umol/L"), cr.mmol.L := data2 * 0.001]
quantile(df$cr.mmol.L[df$Specimen.Unit.Of.Measure %in% c("umol/L")])

# remove abnormally small values; see http://www.mps.com.au/media/1648771/urine%20insight%20final%20low%20res.pdf, although in a different context
# Check the operators
test <- df[cr.mmol.L < 1.8]
table(test$Operator, exclude = NULL)
# only < & =, so remove all
df <- df[cr.mmol.L >= 1.8]

# plausible values

range(df$cr.mmol.L)
subset(df, cr.mmol.L >= 30000) # everything above 40,000 looks very large, remove
subset(df, patid == 4692147) # cr.mmol.L > 35300; rest much smaller, so seems a typo; remove that one too
df <- df[cr.mmol.L <= 30000]

# final clean & save

df_cr <- .clean_df_post(df, colname = "cr.mmol.L", mean_per_day = FALSE)

### combine albumin & creatinine

df <- merge(df_alb, df_cr, by = c("patid", "eventdate"), suffix = c("_alb", "_cr"))
df[, val := alb.mg.L / cr.mmol.L]
df[, alb.mg.L := NULL]
df[, cr.mmol.L := NULL]

# interpretable values

df <- df[(Operator_alb %in% c("=", "~") & Operator_cr %in% c("=", "~")) |
           (Operator_alb %in% c("<", "<=", "=", "~")) & (Operator_cr %in% c(">", ">=", "=", "~") & val < lb) |
           (Operator_alb %in% c("<")) & (Operator_cr %in% c(">", ">=", "=") & val == lb) |
           (Operator_alb %in% c("<", "<=", "=", "~")) & (Operator_cr %in% c(">") & val == lb) |
           (Operator_alb %in% c(">", ">=", "=", "~")) & (Operator_cr %in% c("<", "<=", "=", "~") & val > ub) |
           (Operator_alb %in% c(">")) & (Operator_cr %in% c("<", "<=", "=", "~") & val == ub) |
           (Operator_alb %in% c(">", ">=", "=", "~")) & (Operator_cr %in% c("<") & val == ub) |
           val == 0]

# plausible values

subset(df, val >= 10000) # nothing obviously inconsistent; leave as is
df[, Operator_alb := NULL]
df[, Operator_cr := NULL]

# assign stage 

df <- .assign_stage_cts(df, lb, ub)
table(df$stage, exclude = NULL)

# add to the output dataset

df[, test := "acr.mg.mmol"]
df[, lab := "acr, ratio"]

output <- rbind(output, df)
output

# check
output[, .(mean = mean(val)), by = "lab"]
output[, .(prob25 = quantile(val, probs = 0.25)), by = "lab"]
output[, .(prob75 = quantile(val, probs = 0.75)), by = "lab"]
output[, .(max = max(val)), by = "lab"]
output[, .(min = min(val)), by = "lab"]

# seems fairly consistent, so assume correct

######################################################################
### Numeric values: Protein-to-creatinine ratio
######################################################################

### Scenario 1: stand-alone PCR test ###

# initial cleaning
df <- .clean_df(filename = "df_pcr.Rdata")

### units & conversion

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# check data2 = 0
test <- df[data2 == 0]
table(test$Specimen.Unit.Of.Measure)

# keep only units vaguely similar to mg/mmol or mg/g
# KEEP units that have (creat) in them
df <- df[(Specimen.Unit.Of.Measure %in% c("g/mol", "mg/mg", "mg/mmol", "mg/mmol(creat)")) |
           (data2 == 0 & Specimen.Unit.Of.Measure %in% c("g/mol", "mg/mmol", "mg/mmol(creat)", "No Data Entered", "ratio"))]

# convert all to mg/mmol
# in fact, g/mol is the same as mg/mmol
df[, pcr.mg.mmol := NaN]

df[Specimen.Unit.Of.Measure %in% c("g/mol", "mg/mmol", "mg/mmol(creat)"), pcr.mg.mmol := data2]
quantile(df$pcr.mg.mmol[df$Specimen.Unit.Of.Measure %in% c("g/mol")])
quantile(df$pcr.mg.mmol[df$Specimen.Unit.Of.Measure %in% c("mg/mmol")]) 
quantile(df$pcr.mg.mmol[df$Specimen.Unit.Of.Measure %in% c("mg/mmol(creat)")])
# possibly some values abnormally high (or low!) but nothing consistent

df[Specimen.Unit.Of.Measure %in% c("mg/mg"), pcr.mg.mmol := data2 * 1000 * 0.113]
quantile(df$pcr.mg.mmol[df$Specimen.Unit.Of.Measure %in% c("mg/mg")])

df[data2 == 0, pcr.mg.mmol := 0]

### interpretable values

table(df$Operator, exclude = NULL)
lb <- 15
ub <- 50
df <- .filter_lb_up(df, colname = "pcr.mg.mmol", lb = lb, ub = ub)
table(df$Operator, exclude = NULL)

### plausible values

range(df$pcr.mg.mmol)
subset(df, pcr.mg.mmol >= 20000) # the value of >40000 abnormally high
df <- df[pcr.mg.mmol <= 40000]
df[, Operator := NULL]

### final clean

df <- .clean_df_post(df = df, colname = "pcr.mg.mmol", mean_per_day = FALSE)
df <- rename(df, val = pcr.mg.mmol)

# assign stage 

df <- .assign_stage_cts(df, lb, ub)
table(df$stage, exclude = NULL)

# add to the output

df[, test := "pcr.mg.mmol"]
df[, lab := "pcr, standalone"]
output <- rbind(output, df)

### Scenario 2: separate urine protein & urine creatinine tests ###

### urine protein concentrations

df <- .clean_df(filename = "df_prot.Rdata")

# units & conversion
table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# check data2 = 0
test <- df[data2 == 0]
table(test$Specimen.Unit.Of.Measure)

# leave only units vaguely look like mg/L
df <- df[(Specimen.Unit.Of.Measure %in% c("g/dL", "g/L", "mg/dL", "mg/L")) |
           (data2 == 0 & Specimen.Unit.Of.Measure %in% c("g/L", "mg/L", "No Data Entered"))]

# convert into mg/L

df[, prot.mg.L := NaN]

df[Specimen.Unit.Of.Measure %in% c("g/dL"), prot.mg.L := data2 * 10 * 1000]
quantile(df$prot.mg.L[df$Specimen.Unit.Of.Measure %in% c("g/dL")]) # consistently too high; remove
df <- df[Specimen.Unit.Of.Measure != "g/dL"]

df[Specimen.Unit.Of.Measure %in% c("g/L"), prot.mg.L := data2 * 1000]
quantile(df$prot.mg.L[df$Specimen.Unit.Of.Measure %in% c("g/L")]) # some values too high, but some could be plausible, so leave for now & filter later

df[Specimen.Unit.Of.Measure %in% c("mg/dL"), prot.mg.L := data2 * 10]
quantile(df$prot.mg.L[df$Specimen.Unit.Of.Measure %in% c("mg/dL")]) # consistently too high; remove
df <- df[Specimen.Unit.Of.Measure != "mg/dL"]

df[Specimen.Unit.Of.Measure %in% c("mg/L"), prot.mg.L := data2]
quantile(df$prot.mg.L[df$Specimen.Unit.Of.Measure %in% c("mg/L")])

df[data2 == 0, prot.mg.L := 0]

# plausible values

range(df$prot.mg.L)
subset(df, prot.mg.L >= 100000) # one value abnormally high
subset(df, patid %in% subset(df, prot.mg.L >= 100000)$patid) # all large values seem a misprint - or rather wrong units; remove
df <- df[prot.mg.L < 100000]

# save
df_prot <- .clean_df_post(df, colname = "prot.mg.L", mean_per_day = FALSE)

### urine creatinine concentrations

# already done at the albuminuria stage

### combine protein & creatinine

df <- merge(df_prot, df_cr, by = c("patid", "eventdate"), suffix = c("_prot", "_cr"))
df[, val := prot.mg.L / cr.mmol.L]
df[, prot.mg.L := NULL]
df[, cr.mmol.L := NULL]

# interpretable values

df <- df[(Operator_prot %in% c("=", "~") & Operator_cr %in% c("=", "~")) |
           (Operator_prot %in% c("<", "<=", "=", "~")) & (Operator_cr %in% c(">", ">=", "=", "~") & val < lb) |
           (Operator_prot %in% c("<")) & (Operator_cr %in% c(">", ">=", "=") & val == lb) |
           (Operator_prot %in% c("<", "<=", "=", "~")) & (Operator_cr %in% c(">") & val == lb) |
           (Operator_prot %in% c(">", ">=", "=", "~")) & (Operator_cr %in% c("<", "<=", "=", "~") & val > ub) |
           (Operator_prot %in% c(">")) & (Operator_cr %in% c("<", "<=", "=", "~") & val == ub) |
           (Operator_prot %in% c(">", ">=", "=", "~")) & (Operator_cr %in% c("<") & val == ub) |
           val == 0]

### plausible values

range(df$val)
subset(df, val >= 3000)
subset(df, patid %in% subset(df, val >= 5000)$patid) # all seems consistent; one patient with only one test of 6300, but then there are people iwth similar values
# leave all as is
df[, Operator_prot := NULL]
df[, Operator_cr := NULL]

# assign stage 

df <- .assign_stage_cts(df, lb, ub)
table(df$stage, exclude = NULL)

# add to the output dataset

df[, test := "pcr.mg.mmol"]
df[, lab := "pcr, ratio"]
output <- rbind(output, df)

### check consistency between the two PCR definitions

alpha <- output[lab %in% c("pcr, standalone", "pcr, ratio")]
alpha[, .(mean = mean(val)), by = "lab"]
alpha[, .(prob25 = quantile(val, probs = 0.25)), by = "lab"]
alpha[, .(prob75 = quantile(val, probs = 0.75)), by = "lab"]
alpha[, .(max = max(val)), by = "lab"]
alpha[, .(min = min(val)), by = "lab"]

# standalone pcr slightly higher, but otherwise no massive difference; assume correct

######################################################################
### Numeic values: Albumin excretion rate
######################################################################

### Scenario 1: stand-alone AER test ###

# initial cleaning
df <- .clean_df(filename = "df_aer.Rdata")

# units & conversion

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# check data2 = 0
test <- df[data2 == 0]
table(test$Specimen.Unit.Of.Measure)

# remove only units vaguely similar to mg/24hrs
df <- df[(Specimen.Unit.Of.Measure %in% c("g/d", "ug/min")) |
           (data2 == 0 & Specimen.Unit.Of.Measure %in% c("No Data Entered"))]
#nrow(df[data2 == 0]) - nrow(test)

#convert all to mg/day
df[, aer.mg.day := NaN]

df[Specimen.Unit.Of.Measure %in% c("g/d"), aer.mg.day := data2 * 1000]
quantile(df$aer.mg.day[df$Specimen.Unit.Of.Measure %in% c("g/d")])

df[Specimen.Unit.Of.Measure %in% c("ug/min"), aer.mg.day := data2 * 1.44]
quantile(df$aer.mg.day[df$Specimen.Unit.Of.Measure %in% c("ug/min")]) 

df[data2 == 0, aer.mg.day := 0]

# interpretable values

table(df$Operator, exclude = NULL)
lb <- 30
ub <- 300
df <- .filter_lb_up(df, colname = "aer.mg.day", lb = lb, ub = ub)

# plausible values

range(df$aer.mg.day)
subset(df, aer.mg.day >= 1000)
subset(df, patid %in% subset(df, aer.mg.day >= 1000)$patid) # only one test per patient; however, nothing too abnormal
# leave as is

# final clean

df <- .clean_df_post(df = df, colname = "aer.mg.day", mean_per_day = FALSE)
df[, Operator := NULL]
df <- rename(df, val = aer.mg.day)

# assign stage 

df <- .assign_stage_cts(df, lb, ub)
table(df$stage, exclude = NULL)

# add to the output

df[, test := "aer.mg.day"]
df[, lab := "aer, standalone"]
output <- rbind(output, df)

######################################################################
### Numeric values: Protein excretion rate
######################################################################

### Scenario 1: stand-alone PER test ###

# initial cleaning
df <- .clean_df(filename = "df_per.Rdata")

# units & conversion

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# check data2 = 0
test <- df[data2 == 0]
table(test$Specimen.Unit.Of.Measure)
# remove only units vaguely similar to mg/24hrs
df <- df[(Specimen.Unit.Of.Measure %in% c("g/d", "ug/d")) |
           (data2 == 0 & Specimen.Unit.Of.Measure %in% c("g/d", "No Data Entered"))]

# convert all to mg/day
df[, per.mg.day := NaN]

df[Specimen.Unit.Of.Measure %in% c("g/d"), per.mg.day := data2 * 1000]
quantile(df$per.mg.day[df$Specimen.Unit.Of.Measure %in% c("g/d")])

df[Specimen.Unit.Of.Measure %in% c("ug/d"), per.mg.day := data2 * 0.001]
quantile(df$per.mg.day[df$Specimen.Unit.Of.Measure %in% c("ug/d")]) # very low; remove
df <- df[Specimen.Unit.Of.Measure != "ug/d"]

df[data2 == 0, per.mg.day := 0]

# interpretable values

table(df$Operator, exclude = NULL)
lb <- 150
ub <- 500
df <- .filter_lb_up(df, colname = "per.mg.day", lb = lb, ub = ub)

# plausible values

subset(df, per.mg.day >= 9000) # nothing too much out of order; leave as is
subset(df, patid %in% subset(df, per.mg.day >= 9000)$patid) # some values possibly too high, but others look plausible; leave as is

# final clean
df <- .clean_df_post(df = df, colname = "per.mg.day", mean_per_day = FALSE)
df <- rename(df, val = per.mg.day)

# assign stage 

df <- .assign_stage_cts(df, lb, ub)
table(df$stage, exclude = NULL)

# add to the output

df[, test := "per.mg.day"]
df[, lab := "per, standalone"]
output <- rbind(output, df)

######################################################################
### Numerical: check combined results
######################################################################

output[, .(mean = mean(val)), by = "lab"]
output[, .(prob25 = quantile(val, probs = 0.25)), by = "lab"]
output[, .(prob75 = quantile(val, probs = 0.75)), by = "lab"]
output[, .(max = max(val)), by = "lab"]
output[, .(min = min(val)), by = "lab"]

table(output$stage, output$test, exclude = NULL)

### save

output_test_num <- output
setwd(outputDir)
save(output_test_num, file = "albuminuria_test_num.Rdata")
load("albuminuria_test_num.Rdata")
head(output_test_num)
table(output_test_num$stage, output_test_num$test)


######################################################################
### Categorical tests
######################################################################

# load test data

setwd(sourceDir_R_raw)
df <- get(load("df_alb_test.Rdata"))

# add readterms

setwd(sourceDir_codelists)
codes <- fread("albuminuria_all_BF_final.csv")
codes <- codes[, .(medcode, readterm)]
df <- merge(df, codes, by = "medcode", all.x = TRUE)
df[, medcode := NULL]

# decipher field values

df <- merge(df, entity[, .(enttype, data1, data1_lkup, data2, data2_lkup, 
                           data3, data3_lkup, data4, data4_lkup)], by = "enttype", suffix = c("", "_def"), all.x = TRUE)
table(df$data1_def)

### numerical output ###

alpha <- df[data1_def == "Operator"]
sort(unique(alpha$readterm))

# check the categorical diagnoses that we haven't picked up?
beta <- alpha[readterm %in% c("[D]Microalbuminuria", "[D]Proteinuria", "Isolated proteinuria with specified morphological lesion", 
              "Proteinuria")]
table(beta$readterm, beta$enttype, exclude = NULL)

# Microalbuminuria: 304 values against "Urine microalbumin"; this was assumed not to be the ratio & extracted in the "albumin" codelist, so not interpretable
# 3 proteinuria codes against 884 entries of "Other laboratory tests"; presumably pcr, but assume not interpretable

### categorical output ###

alpha <- df[data1_def == "Qualifier"]

# only leave meaningful columns

alpha <- alpha[, .(patid, eventdate, readterm, enttype, data1)]

# add test value information

alpha <- rename(alpha, Code = data1)
alpha <- merge(alpha, TQU, by = "Code", all.x = TRUE)
alpha[, Code := NULL]

### remove values that cannot be interpreted or seem wrongly coded

# test qualifiers that cannot be interpreted
table(alpha$Test.Qualifier, exclude = NULL)
alpha <- alpha[!(Test.Qualifier %in% c("High", "Low", "Not examined", "Outside ref range", "Potential Abnormal", "Potentially abnormal"))]

# readterms with no value
table(alpha$readterm, exclude = NULL)
alpha <- alpha[!(readterm %in% c("Urine creatinine", "Urine microalbumin positive", "Urine protein test not done"))]

# wrong enttypes
table(alpha$enttype, exclude = NULL)
alpha <- alpha[!(enttype %in% c(229, 263, 264))]
# 229: pregnancy test
# 263: ribs x-ray
# 264: shoulders x-ray
nrow(alpha)

# TODO: additional check for 287: protein & 431: urine dipstick for protein?
# leave as is for now & use readcodes & test values instead

# remove entries with Data Not Entered unless the readterm is self-explanatory
alpha <- rename(alpha, val = Test.Qualifier)
readterms <- c("Urine microalbumin negative",
               "Urine protein normal", 
               "Urine protein test = +", "Urine protein test = ++", "Urine protein test = +++", "Urine protein test = ++++", 
               "Urine protein test = trace",
               "Urine protein test negative")
alpha <- alpha[(val != "Data Not Entered") |
                (val == "Data Not Entered" & readterm %in% readterms)]
table(alpha$val)

### create intepretable combinations

# possible proteinuria readterm & result of a urine protein test
tests <- expand.grid(c("[D]Proteinuria", "[D]Proteinuria NOS", 
                       "24 hour urine protein excretion test", "24 hour urine protein output",
                       "Proteinuria", "Random urine protein level",
                       "Urine dipstick for protein", "Urine protein",
                       "Urine protein abnormal", "Urine protein level",
                       "Urine protein NOS", "Urine protein test", 
                       "Urine protein test NOS", "Urine total protein"),
                     c("+", "++", "+++", "Negative", "Nil", "Normal", "Trace"))
# defnite proteinuria readterm & positive result of a urine protein test
tests <- rbind(tests, 
               expand.grid(c("Type 2 diabetes mellitus with persistent proteinuria"),
                           c("+", "++", "+++", "Trace")))
# suggestive readterms & coresponding test result
tests <- rbind(tests,
            expand.grid(c("Urine microalbumin negative"),
                        c("Negative", "Nil", "Normal", "Data Not Entered")),
            expand.grid(c("Urine protein abnormal"),
                        c("+", "++", "+++", "Trace")),
            expand.grid(c("Urine protein normal"),
                        c("Negative", "Nil", "Normal", "Data Not Entered")),
            expand.grid(c("Urine protein test = +"), 
                        c("Abnormal", "Positive", "+", "Data Not Entered")),
            expand.grid(c("Urine protein test = ++"), 
                        c("Abnormal", "Positive", "++", "Data Not Entered")),
            expand.grid(c("Urine protein test = +++"), 
                        c("Abnormal", "Positive", "+++", "Data Not Entered")),
            # NB: no "++++" in the val column; the highest is +++ so permit this
            expand.grid(c("Urine protein test = ++++"), 
                        c("Abnormal", "Positive", "+++", "Data Not Entered")),
            expand.grid(c("Urine protein test = trace"),
                        c("Abnormal", "Positive", "Trace", "Data Not Entered")),
            expand.grid(c("Urine protein test negative"), 
                        c("Negative", "Nil", "Normal", "Data Not Entered"))
            )
# non-protein tests & interpretable (= "no") result of a urine protein test
tests <- rbind(tests,
               expand.grid(c("[D]Albuminuria", "Albumin / creatinine ratio", 
                             "Microalbumin excretion rate", "Overnight albumin excretion rate",
                             "Urine albumin", "Urine albumin:creatinine ratio", 
                             "Urine microalbumin", "Urine microalbumin:creatinine ratio"), 
                           c("Negative", "Nil", "Normal")))

tests <- data.table(tests)
tests <- tests[order(Var1)]
tests <- rename(tests, readterm = Var1, val = Var2)

# save
setwd("K:\\SHARP_EE\\Monitoring CKD CHF\\Exploratory analyses\\output\\")
write.csv(tests, file = "cat_tests_scenarios.csv")

# filter out permissible combinations
beta <- merge(alpha, tests, by = c("readterm", "val"))

# checks
gamma <- alpha[patid %in% setdiff(alpha$patid, beta$patid)]
table(gamma$readterm, gamma$val)
write.csv(table(gamma$readterm, gamma$val), file = "cat_test_exclusions.csv")

### assign stage

table(beta$val, exclude = NULL)
beta[, stage := NaN]

beta[val %in% c("Negative", "Nil", "Normal", "Trace"), stage := 1]
beta[val %in% c("+"), stage := 2]
beta[val %in% c("++", "+++"), stage := 3]

gamma <- beta[is.na(stage)]
sort(unique(as.character(gamma$readterm)))

beta[is.na(stage) & readterm %in% c("Urine microalbumin negative", "Urine protein normal", "Urine protein test negative",
                                    "Urine protein test = trace"), stage := 1] 
beta[is.na(stage) & readterm %in% c("Urine protein test = +"), stage := 2]
beta[is.na(stage) & readterm %in% c("Urine protein test = ++", 
                                    "Urine protein test = +++", "Urine protein test = ++++"), stage := 3]

table(beta$stage, exclude = NULL) # all covered

### clean & save

output_test_cat <- beta
setwd(outputDir)
save(output_test_cat, file = "albuminuria_test_cat.Rdata")

######################################################################
### Clinical
######################################################################

# load test data

setwd(sourceDir_R_raw)
df <- get(load("df_alb_clin.Rdata"))

# add readterms

setwd(sourceDir_codelists)
codes <- fread("albuminuria_all_BF_final.csv")
codes <- codes[, .(medcode, readterm)]
df <- merge(df, codes, by = "medcode", all.x = TRUE)
df[, medcode := NULL]

# leave only meaningful fields

alpha <- df[!is.na(data1)]
table(alpha$readterm)
table(alpha$enttype)

# 14: height
# 26: current diabetes status
# 65: diabetic consultation
# 137: allergy and intolerance
# 149: cause of death
# nothing relevant, so discard all this information; only use readterms

df <- df[, .(patid, eventdate, readterm)]

### leave only codes that can be interpreted

sort(unique(df$readterm))

readterms_stage1 <- c("Chronic kidney disease stage 1 without proteinuria",
                      "Chronic kidney disease stage 2 without proteinuria",
                      "Chronic kidney disease stage 3 without proteinuria",
                      "Chronic kidney disease stage 3A without proteinuria",
                      "Chronic kidney disease stage 3B without proteinuria",
                      "Chronic kidney disease stage 4 without proteinuria",
                      "Chronic kidney disease stage 5 without proteinuria",
                      "CKD stage 2 without proteinuria",
                      "CKD stage 3 without proteinuria",
                      "CKD stage 3A without proteinuria",
                      "CKD stage 3B without proteinuria",
                      "CKD stage 4 without proteinuria",
                      "CKD stage 5 without proteinuria",
                      "Urine protein normal",
                      "Urine protein test negative",
                      "Urine microalbumin negative",
                      "Urine protein test = trace"
                      )

readterms_stage2 <- c("Urine protein test = +")

readterms_stage3 <- c("Urine protein test = ++",
                      "Urine protein test = +++", "Urine protein test = ++++")

df <- df[readterm %in% c(readterms_stage1, readterms_stage2, readterms_stage3)]

### assign stage

df[, stage := NaN]
df[readterm %in% readterms_stage1, stage := 1]
df[readterm %in% readterms_stage2, stage := 2]
df[readterm %in% readterms_stage3, stage := 3]
table(df$stage, exclude = NULL)

### save

output_clin <- df
setwd(outputDir)
save(output_clin, file = "albuminuria_clin.Rdata")

######################################################################
### Combine results
######################################################################

# combine thre three datasets

alpha_num <- output_test_num[, .(patid, eventdate, stage)]
alpha_cat <- output_test_cat[, .(patid, eventdate, stage)]
alpha_clin <- output_clin[, .(patid, eventdate, stage)]

df <- rbind(alpha_num, alpha_cat, alpha_clin)
# remove duplicates; although technically they could come from re-testing
df <- distinct(df)

# clean

df <- df[order(patid, eventdate)]
table(df$stage, exclude = NULL)
length(unique(df$patid))
length(unique(alpha_num$patid))

# save

setwd(outputDir)
save(df, file = "df_albuminuria_all.Rdata")

######################################################################
### Calculate yearly medians
######################################################################

### extract suty entry dates

setwd(sourceDir_R)
get(load("2017-10-26 df_base.Rdata"))
df_studyentry <- df_base[, .(patid, studyentry)]

### numerical tests only ###

df_alb <- output_test_num[, .(patid, eventdate, stage)]

# clean
df_alb <- distinct(df_alb)
df_alb <- df_alb[order(patid, eventdate)]
table(df_alb$stage, exclude = NULL)

# put tests in years depending on the study entry
df_alb <- merge(df_alb, df_studyentry, by = "patid") # this also removes patients that are not in df_base
df_alb[, eventyear := ceiling(as.numeric(((eventdate - 3 * 30) - studyentry) / 365.25))]

# For medians take ceiling so that 1.5 is mapped into 2 & 2.5 is mapped into 3 (because it indicates at least one test of a higher stage)
df_alb_peryear <- df_alb[, .(stage = ceiling(median(stage))), by = c("patid", "eventyear")]
df_alb_peryear <- df_alb_peryear[order(patid)]
table(df_alb_peryear$stage, exclude = NULL)

# save
setwd(outputDir)
save(df_alb_peryear, file = "df_alb_peryear_TESTNUM.Rdata")

### categorical tests & diagnoses only ###

df_alb <- rbind(output_test_cat[, .(patid, eventdate, stage)],
                output_clin[, .(patid, eventdate, stage)])

# clean
df_alb <- distinct(df_alb)
df_alb <- df_alb[order(patid, eventdate)]
table(df_alb$stage, exclude = NULL)

# put tests in years depending on the study entry
df_alb <- merge(df_alb,df_studyentry, by = "patid") # this also removes patients that are not in df_base
df_alb[, eventyear := ceiling(as.numeric(((eventdate - 3 * 30) - studyentry) / 365.25))]

# For medians take ceiling so that 1.5 is mapped into 2 & 2.5 is mapped into 3 (because it indicates at least one test of a higher stage)
df_alb_peryear <- df_alb[, .(stage = ceiling(median(stage))), by = c("patid", "eventyear")]
df_alb_peryear <- df_alb_peryear[order(patid)]
table(df_alb_peryear$stage, exclude = NULL)

# save
setwd(outputDir)
save(df_alb_peryear, file = "df_alb_peryear_TESTCAT_CLIN.Rdata")

######################################################################
### Calculate yearly medians (all)
######################################################################
#
#df_alb <- df
#
## extract year of the test
#
#df_alb[, eventyear := as.numeric(format(eventdate, "%Y"))]
#
## For medians take ceiling so that 1.5 is mapped into 2 & 2.5 is mapped into 3 (because it indicates at least one test of a higher stage)
#df_alb_peryear <- df_alb[, .(stage = ceiling(median(stage))), by = c("patid", "eventyear")]
#df_alb_peryear[, eventdate := as.Date(paste0(eventyear, "-07-01"))]
#df_alb_peryear <- df_alb_peryear[, .(patid, eventdate, stage)]
#df_alb_peryear <- df_alb_peryear[order(patid)]
#table(df_alb_peryear$stage, exclude = NULL)
#
## save
#setwd(outputDir)
#save(df_alb_peryear, file = "df_alb_peryear.Rdata")