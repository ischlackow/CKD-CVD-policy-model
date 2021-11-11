### NB To change start date & output filename/directory, go to the second section
### "Define flexible parameters"

######################################################################
### Load packages and files
######################################################################

rm(list = ls())

library(readstata13)
library(dplyr)
library(gdata)
library(data.table)

sourceDir_codelist <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\Codelists\\"
sourceDir_lookup <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Lookups\\TXTFILES\\"
sourceDir_R_raw <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Raw CPRD data\\"
sourceDir_R <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"

### load functions

source("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\code\\functions.R")

### lookups

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

# entity
setwd(sourceDir_R)
entity <- data.table(get(load("entity.Rdata")))
setkey(entity, enttype)

### practice locations

setwd(sourceDir_R)
df_practice <- data.table(get(load("df_practice.Rdata")))

### load basic patient information

setwd(sourceDir_R)
df_pat <- data.table(get(load("df_pat.Rdata")))
df_pat <- df_pat[, .(patid, sex, dob, deathdate)]
setkey(df_pat, patid)

##############################################
## Initial cohort of patients 
##############################################

### add censoring dates

setwd(sourceDir_R)
df_cens <- get(load("df_cens.Rdata"))
setkey(df_cens, patid)

### ONS

setwd(sourceDir_R)
df_death <- data.table(get(load("df_death.Rdata")))
df_death <- df_cens[df_death, nomatch = 0, on = "patid"]
df_death <- df_death[deathdate <= end]
setkey(df_death, patid)

######################################################################
### Define flexible parameters
######################################################################

### provisional study entry date ###

### we are not filtering by practice dates at this stage; only the start date
startdate <- as.Date("2005-01-01")
#startdate <- as.Date("2009-01-01")

### output filename ###

outputDir <- "K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\"
outputFilename <- "df_baseline"

######################################################################
### Initialise the dataset: study entry & CKD stage at entry
######################################################################

# load necessary data

setwd(sourceDir_R)

# earliest date the CKD could have been diagnosed
df_ckd <- get(load("df_ckd.Rdata"))

# all eGFR tests
df_egfr <- get(load("df_egfr.Rdata"))
length(unique(df_egfr$patid))

# leave only tests taken on of after the CKD diagnosis date
# CKD diagnosis date needs to be after the start date
df_egfr <- merge(df_egfr, df_ckd, by = "patid")
length(unique(df_egfr$patid))

# also merge with df_cens to leave only patients with start <= 2013-12-31
df_egfr <- merge(df_egfr, df_cens, by = "patid")
df_egfr[, studyentry := pmax(ckd_date, startdate, start)]
length(unique(df_egfr$patid))

# remove patients with studyentry after end
df_egfr <- df_egfr[studyentry <= end]
length(unique(df_egfr$patid))

### remove patients that had a transplant prior to study entry ###

# data on all transplants
setwd(sourceDir_R_raw)
df <- data.table(get(load("df_transplant_SS.Rdata")))
# first transplant
alpha <- df[order(patid, eventdate)]
alpha <- alpha[, head(.SD, 1), by = patid]
nrow(alpha)

# study entry for those patients
id <- intersect(alpha$patid, df_egfr$patid)
beta <- df_egfr[patid %in% id, .(patid, studyentry)]
beta <- unique(beta)
nrow(beta)

# combine
gamma <- merge(alpha, beta, by = "patid")
gamma <- gamma[eventdate <= studyentry]
id <- gamma$patid
df_egfr <- df_egfr[!(patid %in% id)]
length(id)

### remove women pregnant in 0-12 months before the studyentry ###

# data on all pregnancies
setwd(sourceDir_R_raw)
df <- data.table(get(load("df_pregnancy.Rdata")))
df <- df[, .(patid, eventdate)]
# eligible dates 0-12 months from the diagnosis
df[, d1 := eventdate]
df[, d2 := eventdate + 365]
alpha <- df[, .(patid, d1, d2)]
length(unique(alpha$patid))
alpha

# study entry for those patients
id <- intersect(alpha$patid, df_egfr$patid)
beta <- df_egfr[patid %in% id, .(patid, studyentry)]
beta <- unique(beta)

# combine
gamma <- merge(alpha, beta, by = "patid")
gamma <- gamma[studyentry >= d1 & studyentry <= d2]
id <- gamma$patid
length(unique(id))

# add pregnancy flag by now; as some of the patients are male so presumably a typo, but we want to keep those in
df_egfr[, preg := as.numeric(patid %in% id)]

### initialise dataset & calculate CKD stage at study entry ### 

# calculate the latest non-ckd1 stage taken on or before the ckd_date
df_egfr <- df_egfr[ckd != "ckd1" & eventdate <= studyentry]
df_egfr <- .latest(df_egfr)
df_egfr[, ckd := drop.levels(ckd)]

# clean & initialise the main output dataset
df_base <- df_egfr[, .(patid, pracid, start, end, studyentry, egfr, ckd, preg)]
length(unique(df_base$patid))

### aux information to be used throughout

# ids of relevant patients
id_ckd <- df_base$patid
df_pat <- df_pat[patid %in% id_ckd]
df_cens <- df_cens[patid %in% id_ckd]
df_death <- df_death[patid %in% id_ckd]

# study entry for each patient
df_studyentry <- df_base[, .(patid, studyentry)]
setkey(df_studyentry, patid)

######################################################################
### Cardiovascular disease at baseline
######################################################################

### CPRD

df_cv <- .clean_df(filename = "df_cv.Rdata")
df_cv

filename <- "df_cv.Rdata"
eventdate_colname <- "eventdate"
cens <- FALSE
df <- df_cv

alpha <- .split_df(df_cv)

df_cv_pre <- alpha$dfpre
df_cv_post <- alpha$dfpost
id_CV <- unique(df_cv_pre$patid)

### HES

# load
df_hes <- .clean_df(sourceDir = sourceDir_R, filename = "df_hes.Rdata", eventdate_colname = "discharged")
alpha <- .split_df(df_hes, "discharged")
df_hes_pre <- alpha$dfpre

setwd(sourceDir_codelist)

#ICD10 codes
codes_CV_ICD10 <- fread("CVD_ICD10_final.csv")
icd10_CV <- codes_CV_ICD10$ICD10_code

# extract
df_hes_CV <- df_hes_pre
df_hes_CV[, ICD_3 := substr(ICD, 0, 3)]

# ICD10
df_CV <- df_hes_CV[ICD_3 %in% icd10_CV]
id_CV_hes <- unique(df_CV$patid)

### combine and add information into df_base
id <- union(id_CV, id_CV_hes)
df_base[,  cvd := as.numeric(patid %in% id)]
table(df_base$cvd, exclude = NULL)

######################################################################
### Statin prescription
######################################################################

# 28 days
id <- .get_id_pre_Tx("df_statins.Rdata")
df_base[, statin := as.numeric(patid %in% id)]
table(df_base$statin, df_base$cvd, exclude = NULL)

# 84 days
id_84 <- .get_id_pre_Tx("df_statins.Rdata", days = 84)
df_base[, statin_84 := as.numeric(patid %in% id_84)]
table(df_base$statin, df_base$statin_84, exclude = NULL)

######################################################################
### Proteinuria / albuminuria
######################################################################

setwd(sourceDir_R)

### numeric ###

# load and clean
df_alb <- .clean_df(sourceDir = sourceDir_R, filename = "df_alb_peryear_TESTNUM.Rdata")
df_alb <- df_alb[eventyear <= 0]

# first consider stage 2 or 3 only
alpha <- df_alb[stage >= 2]

# latest stage before study entry
alpha <- alpha[eventyear == max(eventyear), .SD, by = "patid"]

# add to df_base
alpha <- rename(alpha, album_TESTNUM = stage)
df_base <- merge(df_base, alpha[, .(patid, album_TESTNUM)], all.x = TRUE, by = "patid")

# now differentiate between A1 & no measurement
beta <- df_alb[stage == 1]
id <- beta$patid
# unique(df_base$album[df_base$patid %in% id]) # to check
df_base[patid %in% id & is.na(album_TESTNUM), album_TESTNUM := 1]
# the rest would be with no measurements then
table(df_base$album_TESTNUM, exclude = NULL)
table(df_base$album_TESTNUM, df_base$ckd, exclude = NULL)

### categorical (including diagnoses) ###

# load and clean
df_alb <- .clean_df(sourceDir = sourceDir_R, filename = "df_alb_peryear_TESTCAT_CLIN.Rdata")
df_alb <- df_alb[eventyear <= 0]

# first consider stage 2 or 3 only
alpha <- df_alb[stage >= 2]

# latest stage before study entry
alpha <- alpha[eventyear == max(eventyear), .SD, by = "patid"]

# add to df_base
alpha <- rename(alpha, album_TESTCAT = stage)
df_base <- merge(df_base, alpha[, .(patid, album_TESTCAT)], all.x = TRUE, by = "patid")

# now differentiate between A1 & no measurement
beta <- df_alb[stage == 1]
id <- beta$patid
# unique(df_base$album[df_base$patid %in% id]) # to check
df_base[patid %in% id & is.na(album_TESTCAT), album_TESTCAT := 1]
# the rest would be with no measurements then
table(df_base$album_TESTCAT, exclude = NULL)
table(df_base$album_TESTCAT, df_base$ckd, exclude = NULL)
table(df_base$album_TESTNUM, df_base$album_TESTCAT, exclude = NULL)

#alpha <- df_base[cvd == 0 & statin == 0]
#table(alpha$ckd, alpha$album_TESTNUM, exclude = NULL) 

######################################################################
### CRD flag
######################################################################

# add to df_base
id <- .get_id_pre_diag("df_crd.Rdata")
df_base[, crd := as.numeric(patid %in% id)]

table(df_base$crd, exclude = NULL)

######################################################################
### Age and gender
######################################################################

### gender ###

df_base <- df_pat[, .(patid, sex)][df_base, on = "patid"]
table(df_base$sex, exclude = TRUE)

### drop pregnant women ###

table(df_base$preg, exclude = NULL)
df_base <- df_base[preg == 0 & sex == "F" | sex == "M"]
table(df_base$preg, exclude = NULL)
df_base[, preg := NULL]

### age ###

df_base <- df_pat[, .(patid, dob)][df_base, on = "patid"]
# round to the nearest integer
df_base[, age := round(as.numeric(difftime(studyentry, dob, units = "days") / 365.25))]
df_base <- df_base[,  dob := NULL]

nrow(df_base) == length(unique(df_base$patid))

######################################################################
### Ethnicity
######################################################################

### CPRD

# load data
setwd(sourceDir_R_raw)
df_ethnicity <- data.table(get(load("df_ethnicity.Rdata")))
df_ethnicity <- df_ethnicity[patid %in% id_ckd]
nrow(df_ethnicity)
# merge with the categories
setwd(sourceDir_codelist)
ethnicity_cats <- fread("ethnicity_final.csv")

df_ethnicity <- ethnicity_cats[, .(medcode, category)][df_ethnicity, on = "medcode"]
nrow(df_ethnicity)

# remove category 11 (refused to declare)
df <- df_ethnicity[category != 11]

# consider entry closest to study entry
df <- .closest(df)

# now pick min category for each date (very random but SOMETHING needs to be done)
df <- df[order(patid, category)]
df <- df[, head(.SD, 1), by = patid]
df <- df[, .(patid, category)]

### HES ###

setwd(sourceDir_R)
df_hes_pat <- data.table(get(load("df_hes_pat.Rdata")))
nrow(df_hes_pat) == length(unique(df_hes_pat$patid)) # 1 entry per person!
id <- setdiff(df_base$patid, df$patid)
table(df_hes_pat$gen_ethnicity)
alpha <- df_hes_pat[patid %in% id]
alpha[, category  := NaN]
alpha[gen_ethnicity == "Bangladesi", category := 4]
alpha[gen_ethnicity == "Bl_Afric", category := 7]
alpha[gen_ethnicity == "Bl_Carib", category := 6]
alpha[gen_ethnicity == "Bl_Other", category := 10]
alpha[gen_ethnicity == "Chinese", category := 8]
alpha[gen_ethnicity == "Indian", category := 2]
alpha[gen_ethnicity == "Mixed", category := 9]
alpha[gen_ethnicity == "Oth_Asian", category := 5]
alpha[gen_ethnicity == "Other", category := 9]
alpha[gen_ethnicity == "Pakistani", category := 3]
#alpha[gen_ethnicity == "Unknown", category := 11] # leave this out; as this is not the same as refused
alpha[gen_ethnicity == "White", category := 1]
table(alpha$category, alpha$gen_ethnicity, exclude = NULL) # all covered

# combine with df
beta <- alpha[!is.nan(category), .(patid, category)]
df <- rbind(df, beta)

# add category 11 (refused)
id <- unique(setdiff(df_ethnicity$patid, df$patid))

# add to df_base
df <- rename(df, ethnicity = category)
df_base <- df[df_base, on = "patid"]
df_base[patid %in% id, ethnicity := 11]
table(df_base$ethnicity, exclude = NULL)

######################################################################
### Practice location and Townsend deprivation index
######################################################################

setwd(sourceDir_R)
df_imd <- data.table(get(load("df_imd.Rdata")))
table(df_imd$imd2010_5)
df_imd <- df_imd[patid %in% id_ckd]
df_imd <- rename(df_imd, imd = imd2010_5)

df_base <- merge(df_base, df_imd[, .(patid, imd)], by = "patid", all.x = TRUE)
table(df_base$imd, exclude = NULL) # 1 is most affluent; 5 is most deprived 
# NB: not all patients are present in df_imd, so some pracid are missing
# in principle it may be possible to recover this information from the patient files, as these are already patients that have been identified as those from the linked practices
# It is irrelevant for this project though

# add practice location
df_base <- merge(df_base, df_practice[, .(pracid, region_txt)], by = "pracid", all.x = TRUE)
df_base <- rename(df_base, region = region_txt)
table(df_base$region, exclude = NULL)
# west/east split
df_base[, west := region %in% c("North West", "South West", "West Midlands")]
df_base[, east := region %in% c("East Midlands", "East of England", "London", "North East", "South Central", "South East Coast", "Yorkshire & The Humber")]
table(df_base$east, df_base$west, exclude = NULL)

######################################################################
### Systolic blood pressure
######################################################################

#enttype_t <- subset(entity, description == "Blood pressure")$enttype # 1; use clinical file for values
# NB another entity: Ambulatory blood pressure; however the format is not numeric

df_bp <- .clean_df(filename = "df_bp.Rdata")
df_bp <- .split_df(df_bp)$dfpre

# check entity types
filter(entity, enttype %in% unique(df_bp$enttype)) # only 1, which is what we want
table(df_bp$medcode)

# Suitable medcodes:
# 1: O/E - blood pressure reading
# 57: O/E - BP reading
# 100: O/E - blood pressure
# 101: O/E - BP reading very low
# 102: 24hr blood pressure monitoring
# 103: Average 24 hour diastolic pressure; <- KEEP FOR DBP
# 158: New patient screen
# 204: Hypertensive disease
# 351: High blood pressure
# 676: [D] Raised blood pressure reading
# 799: Essential hypertension
# 803: Hypertension screen
# 859: O/E - BP reading normal
# 2666: H/O: hypertension
# 3481: O/E - BP borderline raised
# 4444: Hypertension monitoring
# 5020: O/E - BP reading raised
# 5341: O/E - BP reading low
# 5760: [V]Examination of blood pressure
# 6470: advice
# 10632: White coat hypertension; DO NOT USE
# 13186: Hypertension monitoring
# 13187: CHD monitoring
# 13200: New patient health check
# 14448: Ambulatory blood pressure recording
# 14640: Lying blood pressure reading
# 14641: Standing blood pressure reading
# 14642: Standing blood pressure reading
# 14643: O/E - BP reading very high
# 15126: pre-treatment BP reading
# 18150: Coronary heart disease monitoring administration
# 18418: Sitting systolic blood pressure
# 19905: O/E - BP reading: no postural drop
# 20049: Blood pressure monitoring
# 22476: Blood pressure recorded by patient at home
# 22595: O/E - Diastolic BP reading; <- KEEP FOR DBP
# 23006: O/E - CVS examination NOS
# 23312: O/E - Systolic BP reading
# 25553: O/E - BP labile
# 27271: O/E - BP stable
# 27272: O/E - blood pressure reading NOS
# 27273: O/E - blood pressure decreased
# 27274: O/E - BP bordeline low
# 29261: Average home syslotic blood pressure
# 31305: Average home diastolic blood pressure <- KEEP FPR DBP
# 32208: O/E - Arterial pressure index abnormal
# 34618: Average day interval systolic blood pressure
# 34994: O/E - Arterial presure index normal
# 37242: Standing systolic blood pressure
# 37243: Standing diastolic blood pressure <- KEEP FOR DBP
# 38277: Average 24 hour systolic blood pressure
# 38278: Average night interval systolic blood pressure
# 41052: Ambulatory systolic blood pressure
# 41445: Average day interval diastolic blood pressure <- KEEP FOR DBP
# 42280: Sitting diastolic blood pressurel <- KEEP FOR DBP
# 43282: Lying systolic blood pressure
# 43547: Average night interval diastolic blood pressure <- KEEP FOR DBP
# 48008: Ambulatory diastolic blood pressure <- KEEP FOR DBP
# 65990: Lying diastolic blood pressure <- KEEP FOR DBP
# 94807: Self measured blood pressure reading
# 104203: Referral for ambulatory blood pressure monitoring

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
df <- .latest(df_sbp)

### add to df_base ###

df_base <- df[df_base, on = "patid"]
range(df_base$sbp, na.rm = TRUE)
sum(is.na(df_base$sbp)) / nrow(df_base) # about 6%

######################################################################
### BMI
######################################################################

### height ###

# entity_t <- filter(entity, description == "Height") # entity 14

# leave pre & post values
df <- .clean_df(filename = "df_height.Rdata")

# check entity types
filter(entity, enttype %in% unique(df$enttype)) # only 14, which is what we want
table(df$medcode)
# unsuitable codes:
# 2: O/E - weight
# 6: Influenza vaccination
# 41: Contraception
# 43: Urine protein test
# 71: International normalised ratio
# 78: Asthma
# 81: Asthma monitoring
# 93: Cigarette smoking
# 128: Menorrhagia
# 154: Low back pain
# 204: Hypertensive disease
# 273: Hypothyroidism
# 312: Acute bronchitis
# 415: Acquired toe deformity NOS
# 504: Transient cerebral ischaemia
# 516: Shingles
# 554: Knee joint pain
# 630: Otitis externa NOS
# 636: Anxiey states
# 637: Hyperlipidaemia NOS
# 676: [D]Raised blood pressure reading
# 711: Diabetes mellitus
# 758: Type 2 diabetes mellitus
# 799: Essential hypertension
# 844: Rheumatoid arthritis
# 1001: Chornic obstructive pulmonary disease
# 1268: Paroxysmal atrial fibrillation
# 1430: Angina pectoris
# 1446: Acute exacerbation of chronic obstructive airways disease
# 2074: Glaucoma
# 2122: Pill check
# 2195: Bronchiectasis
# 2281: Acid reflux
# 2542: Whiplash injury
# 2627: Prostatism
# 2666: H/O: hypertension
# 3290: Acquired hypothyroidism
# 4444: Hypertension monitoring
# 5710: Chronic obstructive airways disease NOS
# 6154: O/E - pulse rate
# 7069: Background diabetic retinopathy
# 7378: Asthma management plan given
# 7579: Suspected UTI
# 7622: Smoking cessation advice
# 7735: Cramp
# 8484: Inhaler technique - good
# 9177: Airways obstruction reversible
# 9313: Menopause monitoring
# 10031: Wound care
# 10043: Asthma annual review
# 10674: Influenza vaccination declined
# 10711: Inhaler technique observed
# 11287: Chronic obstructive pulmonary disease annual review
# 11427: Influenza vaccination invitation letter sent
# 12943: Cigar smoker
# 12947: Pipe smoker
# 13173: Asthma not disturbing sleep
# 13174: Asthma not limiting activities
# 13176: Asthma follow-up
# 13187: CHD monitoring
# 14704: Thyroid defficiency
# 16565: Good hypertension control
# 17652: Constipated
# 18150: Coronary heart disease monitoring administration
# 18466: Child height centiles
# 22942: IUD fitting awaited
# 28554: Angina pectoris NOS
# 35655: Influenza vacc.administrat.NOS

df <- df[medcode %in% c(3, 158, 8105, 12362, 15546, 26474, 30763, 35476, 41045, 47516)]

# keep useful information only
df[, cm := as.numeric(data1)]
df[, data1 := NULL]
df <- distinct(df)

# remove implausible values
df <- df[cm >= 0.5 & cm <= 2.30]


# take averages for readings taken on the same day
df <- .dailymean(df, "cm")

# extract entry closest to the study entry
df <- .closest(df)
# take the earliest entry, in case there were two measurements, exactly the same number of days either side of study entry
df <- df[order(patid, eventdate)]
df <- df[, head(.SD, 1), by = patid]
df <- df[, .(patid, cm)]

# add to df_base
df_base <- df[df_base, on = "patid"]
sum(is.na(df_base$cm)) / nrow(df_base) 

### weight ###

# entity_t <- filter(entity, description == "Weight") # entity 13
# format: data1 = weight in kilos; data3 = BMI

# extract data before study entry
df <- .clean_df(filename = "df_weight.Rdata")
df <- .split_df(df)$dfpre

# medcodes
table(df$medcode)

# suitable codes
# 2: O/E: weight
# 126: O/E - Underweight
# 158: new patient screen
# 430: obesity
# 1018: Weight increasing
# 1581: weight symptom
# 2839 O/E - overweight
# 6545: Health education
# 6713: Patient advised to lose weight
# 7984: O/E - obese
# 8041: Weight monitoring
# 8105 Body Mass Index
# 13200: new patient health check
# 13278: Body mass index 30+ - obesity
# 16404: O/E - weight 10-20% over ideal
# 21520: O/E weight NOS
# 22556: Body  mass index 40+ - severely obese
# 23376: O/E - weight within 10% ideal
# 25951: Refer to weight management programme
# 26473: O/E - weight > 20% below ideal
# 28937: Body mass index high k/m2
# 28946: Body mass index normal k/m2
# 29029: O/E weight 10-20% below ideal
# 29538: Follow-up obesity assessment
# 32914: Body mass index low k/m2
# 32974: O/E - weight > 20% over ideal
# 37937: Weight loss from baseline weight
# 38632: treatment of obesity started
# 42309: baseline weight
# 101043: advice given about weight management
# 103499: overweight
# 104002: percentage weight loss

# unsuitable codes
# 1: O/E - blood pressure reading
# 3: O/E - height
# 6: Influenza vaccination
# 27: Alcohol consumption
# 33: Never smoked tobacco
# 41: Contraception

df <- df[medcode %in% c(2, 126, 158, 430, 1018, 1581, 2839, 6545, 6713, 7984, 8041, 8105, 
                        13200, 13278, 16404, 21520, 22556, 23376, 25951, 26473, 28937, 28946, 29029, 29538, 
                        32914, 32974, 37937, 38632, 42309, 101043, 103499, 104002)]
# clean
df[, c("kg", "bmi") := list(as.numeric(data1), as.numeric(data3))]
df <- df[, .(patid, studyentry, eventdate, kg, bmi)]
# remove implausible values
head(df)
range(df$kg)
df[kg < 25, kg := NA]
df[bmi < 10 | bmi > 300, bmi := NA]

# extract latest information

# BMI
df_bmi <- df[!is.na(bmi)]
df_bmi <- .latest(df_bmi)
df_bmi <- .dailymean(df_bmi, "bmi")

# weight
df_kg <- df[!is.na(kg)]
df_kg <- .latest(df_kg)
df_kg <- .dailymean(df_kg, "kg")

# add height
df_cm <- df_base[, .(patid, cm)]

length(setdiff(df_kg$patid, df_bmi$patid)) # almost everyone with weight has bmi
setdiff(df_bmi$patid, df_kg$patid) # almost everyone with BMI has weight

# extract the latest date of these two
df_diff <- merge(df_bmi, df_kg, by = "patid", all = TRUE, suffix = c("_bmi", "_kg"))
nrow(df_diff)
nrow(df_bmi)
df_diff <- df_diff[, .(patid, eventdate_bmi, eventdate_kg)]
df_diff[, diff := as.numeric(eventdate_bmi - eventdate_kg)]
quantile(df_diff$diff, na.rm = TRUE) # mainly 0 but also some positive & negative, which means that the person could have a BMI measurement after weight OR vice versa

### add to df_base ###

# firstly, add those patients who have a bmi measurement, and no later weight measurement OR a bmi measurement only
alpha <- df_diff[eventdate_bmi >= eventdate_kg | is.na(eventdate_kg)]
id <- alpha$patid
beta <- df_bmi[patid %in% id]
df_base <- merge(df_base, beta[, .(patid, bmi)], by = "patid", all.x = TRUE)
sum(is.na(df_base$bmi)) / nrow(df_base)

# now, consider those with a bmi measurement and a later weight measurement OR a kg measurement ONLY
# so a later weight measurement will over-write a previous BMI measurement

alpha <- df_diff[eventdate_bmi < eventdate_kg | is.na(eventdate_bmi)]
id <- alpha$patid
beta <- df_kg[patid %in% id]
beta <- beta[, .(patid, kg)]
beta
df_cm
beta <- merge(beta, df_cm, by = "patid")
beta <- beta[!is.na(cm)]
beta[, bmi := kg / cm^2]
beta <- beta[bmi >= 10 & bmi <= 300]
id <- beta$patid
beta <- beta[, .(patid, bmi)]
beta <- rename(beta, bmi_2 = bmi)

# add to df_base
df_base <- merge(df_base, beta, all.x = TRUE)
df_base[is.na(bmi), bmi := bmi_2]
df_base[, cm := NULL]
df_base[, bmi_2 := NULL]

sum(is.na(df_base$bmi)) / nrow(df_base)
quantile(df_base$bmi, na.rm = TRUE)

nrow(df_base) == length(unique(df_base$patid))

######################################################################
### Smoking
######################################################################

#entity_t <- filter(entity, description == "Smoking") # entity 4
# format: data1 = status; 
# data2 = number of cigarettes per day
# data3 = number of cigars per day
# data4 = ounces of tobacco per day
# data5 = start day
# data6 = stop day
# keep all of these for now

df <- .clean_df(filename = "df_smok.Rdata")
df <- .split_df(df)$dfpre
df <- distinct(df)

# need latest smoking status only

df <- .latest(df)

### number of cigarettes smoked per day

# now add data2, data3 & data4

# convert data2
df[, data2 := as.numeric(data2)]
df[is.na(data2), data2 := 0]
df[data2 >= 300, data2 := 0]
table(df$data2, exclude = NULL)

# convert data3
df[, data3 := as.numeric(data3)]
df[is.na(data3), data3 := 0]
table(df$data3, exclude = NULL)
# NB: Assume 1 cigar = 2 cigarettes. 
# Source: http://www.nhsggc.org.uk/about-us/professional-support-sites/cdm-local-enhanced-services/health-determinants/assess-status/smoking/

# convert data4
df[, data4 := as.numeric(data4)]
df[is.na(df$data4), data4 := 0]
table(df$data4, exclude = NULL)
# NB Assume 1 oz of tobacco = 7 cigarettes per day
# Source: http://www.nhsggc.org.uk/about-us/professional-support-sites/cdm-local-enhanced-services/health-determinants/assess-status/smoking/

# add up
df[, n_smok := data2 + 2 * data3 + 7 * data4]

### determine smoking status

# add codes
table(df$data1, exclude = NULL)
YND <- rename(YND, data1 = Code, status = Smoke.drink.Status)
df[, data1 := as.integer(data1)]
df <- merge(df, YND, by = "data1")
df[, stopdate := as.Date(data6, format = "%d/%m/%Y")]
df <- df[, .(patid, n_smok, status, stopdate)]
df <- df[order(patid)]

df[, category := NaN]

# extract never smokers
beta <- df[n_smok == 0 & is.na(stopdate) & status %in% c("No", "Data Not Entered")]
id <- unique(beta$patid)
df[patid %in% id, category := 1]

# extract ex-smokers
beta <- df[n_smok == 0 & (!is.na(stopdate) | status == "Ex")]
# NB: this picks up those with stopdate & status = "Data Not Entered"
id <- unique(beta$patid)
df[patid %in% id, category := 2]

# current smokers
beta <- df[n_smok > 0 | (n_smok == 0 & is.na(stopdate) & status %in% c("Yes"))]

# light smokers
gamma <- beta[n_smok <= 9]
id <- unique(gamma$patid)
df[patid %in% id, category := 3]

# moderate smokers
gamma <- beta[n_smok >= 10 & n_smok <= 19]
id <- unique(gamma$patid)
df[patid %in% id, category := 4]

# heavy smokers
gamma <- beta[n_smok >= 20]
id <- unique(gamma$patid)
df[patid %in% id, category := 5]

# check
table(df$category, exclude = NULL)

# double check that each patient only has one record
nrow(df) == length(unique(df$patid)) # yes

### add to df_base
df <- df[, .(patid, category)]
df <- rename(df, smoking = category)
df_base <- merge(df_base, df, by = "patid", all.x = TRUE)
table(df_base$smoking, exclude = NULL)

nrow(df_base) == length(unique(df_base$patid))

######################################################################
### Family history of CHD in a first-degree relative < 60
######################################################################

# NB: consider all entries, not just those pre study entry
# All consider all fmaily history unless it says explicitly it's not premature

### records from the READ codes

df <- .clean_df(filename = "df_fh_chd.Rdata")
id <- unique(df$patid)

### additional records from entity 87 & generic READ code of Family history

setwd(sourceDir_codelist)
alpha <- fread("CHD_final.csv")
alpha <- alpha[included == 1]
medcodes <- alpha$medcode

df <- .clean_df(filename = "df_chd_2.Rdata")
df <- df[, .(patid, medcode, data1)]
df[, data1 := as.numeric(data1)]
df <- df[medcode == 1554 & data1 %in% medcodes]
id_2 <- df$patid

### combine & add to df_base

id <- union(id, id_2)
df_base[, fh_chd := as.numeric(patid %in% id)]
table(df_base$fh_chd, exclude = NULL)

######################################################################
### Diabetes
######################################################################

### diagnosis

# all pre-study entry records
df <- .clean_df(filename = "df_diab_diag.Rdata")
df <- .split_df(df)$dfpre

# extract first ever record; assuming this is the diagnosis date
df[, medcode := NULL]
df <- distinct(df)
df <- df[order(patid, eventdate)]
df <- df[, head(.SD, 1), by = patid]
nrow(df) == length(unique(df$patid)) # yes, one record per patient

# calculate age at this date
df <- merge(df, df_pat[, .(patid, dob)], by = "patid")
df[, age := as.numeric(df$eventdate - df$dob)/365.25]

### insulin treatment

# 24 days
id_Tx <- .get_id_pre_Tx("df_ins.Rdata")

# 84 days
id_Tx_84 <- .get_id_pre_Tx("df_ins.Rdata", days = 84)

### combine information and determine type of diabetes

# 28 days
df[, ins := as.numeric(df$patid %in% id_Tx)]
df[, type := ifelse(age <= 35 & df$ins == 1, 1, 2)]
subset(df, age < 0) # age up to -0.5, so this is just the wrong date of birth; leave everything as is

# 84 days
df[, ins_84 := as.numeric(df$patid %in% id_Tx_84)]
df[, type_84 := ifelse(age <= 35 & df$ins_84 == 1, 1, 2)]

# merge with df_base
df <- df[, .(patid, type, type_84)]
df <- rename(df, diab = type, diab_84 = type_84)
df_base <- merge(df_base, df, by = "patid", all.x = TRUE)
table(df_base$diab, exclude = NULL)
df_base[is.na(diab), diab := 0]
table(df_base$diab, exclude = NULL)

table(df_base$diab, df_base$diab_84, exclude = NULL)

#nrow(df_base) == length(unique(df_base$patid))

######################################################################
### Treated hypertension
######################################################################

### diagnosis

id_diag <- .get_id_pre_diag("df_hyp_diag.Rdata")

### prescription of anti-hypertensives

# 28 days
id_Tx <- .get_id_pre_Tx("df_antihyp.Rdata")

# 84 days
id_Tx_84 <- .get_id_pre_Tx("df_antihyp.Rdata", days = 84)

### combine and add information to df_base

id <- intersect(id_diag, id_Tx)
df_base[, tr_hyp := as.numeric(patid %in% id)]
table(df_base$tr_hyp, exclude = NULL)

id_84 <- intersect(id_diag, id_Tx_84)
df_base[, tr_hyp_84 := as.numeric(patid %in% id_84)]

#### On hypertensives at baseline (84 day presecriptions)
df_base[, antihyp := as.numeric(patid %in% id_Tx_84)]

table(df_base$tr_hyp, df_base$tr_hyp_84, exclude = NULL)
table(df_base$tr_hyp_84, df_base$antihyp)

######################################################################
### RA
######################################################################

id <- .get_id_pre_diag("df_ra.Rdata")
df_base[, ra := as.numeric(patid %in% id)]
table(df_base$ra, exclude = NULL)

######################################################################
### Atrial fibrillation
#######################################################################

id <-  .get_id_pre_diag("df_af.Rdata")
df_base[, af := as.numeric(patid %in% id)]
table(df_base$af, exclude = NULL)

######################################################################
### CKD - stage 4 or 5 OR impaired renal function
######################################################################

# stage 4 & 5
id <- .get_id_pre_diag("df_ckd45.Rdata")

# major chronic renal disease
id2 <- .get_id_pre_diag("df_crd.Rdata")

# combine
intersect(id, id2) # very few; let's make it OR
id <- union(id, id2)

df_base[, ckd45_or_crd := as.numeric(patid %in% id)]
table(df_base$ckd45_or_crd, exclude = NULL)

######################################################################
### CKD - stage 3, 4, 5 OR impaired renal function
######################################################################

### extract information on stage 3 & combine with ckd45_or_crd

id3 <- .get_id_pre_diag("df_ckd3.Rdata")
id3A <- .get_id_pre_diag("df_ckd3A.Rdata")
id3B <- .get_id_pre_diag("df_ckd3B.Rdata")
id <- unique(c(id3, id3A, id3B, df_base$patid[df_base$ckd45_or_crd == 1]))

df_base[, ckd345_or_crd := as.numeric(patid %in% id)]
table(df_base$ckd345_or_crd, exclude = NULL)

######################################################################
### CKD stage based on the READ code
######################################################################

# latest CKD 1
df_ckd1 <- .clean_df(filename = "df_ckd1.Rdata")
df_ckd1 <- .split_df(df_ckd1)$dfpre
df_ckd1[, stage := 1]
df_ckd1 <- df_ckd1[, .(patid, eventdate, stage)]
df_ckd1 <- .latest(df_ckd1)

# latest CKD 2
df_ckd2 <- .clean_df(filename = "df_ckd2.Rdata")
df_ckd2 <- .split_df(df_ckd2)$dfpre
df_ckd2[, stage := 2]
df_ckd2 <- df_ckd2[, .(patid, eventdate, stage)]
df_ckd2 <- .latest(df_ckd2)

# latest CKD 3A
df_ckd3A <- .clean_df(filename = "df_ckd3A.Rdata")
df_ckd3A <- .split_df(df_ckd3A)$dfpre
df_ckd3A[, stage := "3A"]
df_ckd3A <- df_ckd3A[, .(patid, eventdate, stage)]
df_ckd3A <- .latest(df_ckd3A)

# CKD 3B
df_ckd3B <- .clean_df(filename = "df_ckd3B.Rdata")
df_ckd3B <- .split_df(df_ckd3B)$dfpre
df_ckd3B[, stage := "3B"]
df_ckd3B <- df_ckd3B[, .(patid, eventdate, stage)]
df_ckd3B <- .latest(df_ckd3B)

# CKD 3
df_ckd3 <- .clean_df(filename = "df_ckd3.Rdata")
df_ckd3 <- .split_df(df_ckd3)$dfpre
df_ckd3[, stage := 3]
df_ckd3 <- df_ckd3[, .(patid, eventdate, stage)]
df_ckd3 <- .latest(df_ckd3)

# CKD 4
df_ckd4 <- .clean_df(filename = "df_ckd4.Rdata")
df_ckd4 <- .split_df(df_ckd4)$dfpre
df_ckd4[, stage := 4]
df_ckd4 <- df_ckd4[, .(patid, eventdate, stage)]
df_ckd4 <- .latest(df_ckd4)

# CKD 5
df_ckd5 <- .clean_df(filename = "df_ckd5.Rdata")
df_ckd5 <- .split_df(df_ckd5)$dfpre
df_ckd5[, stage := 5]
df_ckd5 <- df_ckd5[, .(patid, eventdate, stage)]
df_ckd5 <- .latest(df_ckd5)

# combine & find the latest diagnosis
df <- rbind(df_ckd1, df_ckd2, df_ckd3A, df_ckd3B, df_ckd3, df_ckd4, df_ckd5)
df <- .latest(df)
df[, ckd_READ := factor(stage, levels = c(1, 2, "3A", "3B", 3, 4, 5))]

# add to df_base
df_base <- merge(df_base, df[, .(patid, ckd_READ)], by = "patid", all.x = TRUE)
table(df_base$ckd_READ, df_base$ckd, exclude = NULL)
table(df_base$ckd_READ, df_base$ckd345_or_crd, exclude = NULL)

######################################################################
### Systolic blood pressure variability
######################################################################

### extract data from df_sbp ###

# all BP values in the five years before the study entry
df <- merge(df_sbp, df_studyentry, by = "patid")
df <- df[studyentry - eventdate <= 365.25 * 5]
df <- distinct(df)
nrow(df)

# leave only patients with at least two measurements
alpha <- df[, .N, by = patid]
alpha <- alpha[ N >= 2]
id <- alpha$patid
df <- df[patid %in% id]

# calculate sd
df <- df[, .(sbp_sd = sd(sbp)), by = patid]

### add to df_base ###

df_base <- merge(df_base, df, by = "patid", all.x = TRUE)
quantile(df_base$sbp_sd, na.rm = TRUE)
# replace 0 with 0.1 [checks that this seems okay performed earlier]
# this does indeed seem to correspond to patients whose records all contained the same value
# this could be genuine, a misprint, or perhaps a value carried over from the previous visit
# but we cannot tell so assume correct as is.
df_base[sbp_sd == 0, sbp_sd := 0.1]
quantile(df_base$sbp_sd, na.rm = TRUE)
sum(is.na(df_base$sbp_sd)) / nrow(df_base) # 7% missing
nrow(df_base) == length(unique(df_base$patid))

######################################################################
### Migraine
######################################################################

id <- .get_id_pre_diag("df_migraine.Rdata")
df_base[, migraine := as.numeric(patid %in% id)]
table(df_base$migraine, exclude = NULL)

######################################################################
### Corticosteroid use
######################################################################

# 28 days
id_Tx <- .get_id_pre_Tx("df_steroids.Rdata")
df_base[, steroids := as.numeric(patid %in% id_Tx)]
table(df_base$steroids)

# 84 days
id_Tx_84 <- .get_id_pre_Tx("df_steroids.Rdata", days = 84)
df_base[, steroids_84 := as.numeric(patid %in% id_Tx)]


######################################################################
### Systemic lupus erythematosus
######################################################################

id <- .get_id_pre_diag("df_lupus.Rdata")
df_base[, lupus := as.numeric(patid %in% id)]
table(df_base$lupus, exclude = NULL)

######################################################################
### Second generation "atypical" antipsychotic use
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx("df_antipsych.Rdata")
df_base[, antipsych := as.numeric(patid %in% id_Tx)]
table(df_base$antipsych, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx("df_antipsych.Rdata", days = 84)
df_base[, antipsych_84 := as.numeric(patid %in% id_Tx)]
table(df_base$antipsych, df_base$antipsych_84, exclude = NULL)

######################################################################
### Severe mental illness
######################################################################

id <- .get_id_pre_diag("df_mental.Rdata")
df_base[, mental := as.numeric(patid %in% id)]
table(df_base$mental, exclude = NULL)

######################################################################
### HIV or AIDS
######################################################################

id <- .get_id_pre_diag("df_hiv.Rdata")
df_base[, hiv := as.numeric(patid %in% id)]
table(df_base$hiv, exclude = NULL)

######################################################################
### Erectile dysfunction
######################################################################

# diagnosis
id_diag <- .get_id_pre_diag("df_erectile_diag.Rdata")

# treatment

# 28 days
id_Tx <-  .get_id_pre_Tx("df_erectile_Tx.Rdata")

# combine and add information to df_base
id <- union(id_diag, id_Tx)

df_base[, impotence := as.numeric(patid %in% id)]
# remove for females
df_base[sex == "F", impotence := 0]
table(df_base$impotence, exclude = NULL)
table(df_base$impotence, df_base$sex, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx("df_erectile_Tx.Rdata", days = 84)

# combine and add information to df_base
id_84 <- union(id_diag, id_Tx_84)

df_base[, impotence_84 := as.numeric(patid %in% id_84)]
# remove for females
df_base[sex == "F", impotence_84 := 0]
table(df_base$impotence, df_base$impotence_84, exclude = NULL)



######################################################################
### Diastolic blood pressure
######################################################################

# check entity types
filter(entity, enttype %in% unique(df_bp$enttype)) # only 1, which is what we want
table(df_bp$medcode)
# Medcodes:
# 1: O/E - blood pressure reading
# 57: O/E - BP reading
# 100: O/E - blood pressure
# 101: O/E - BP reading very low
# 102: 24hr blood pressure monitoring
# 103: Average 24 hour diastolic pressure; <- KEEP FOR DBP
# 158: New patient screen
# 351: High blood pressure
# 676: [D] Raised blood pressure reading
# 799: Essential hypertension
# 803: Hypertension screen
# 859: O/E - BP reading normal
# 2379: Seen in diabetic clinic
# 3481: O/E - BP borderline raised
# 4444: Hypertension monitoring
# 5020: O/E - BP reading raised
# 5341: O/E - BP reading low
# 5760: [V]Examination of blood pressure
# 6470: advice
# 7917: Follow-up arranged
# 10632: White coat hypertension; DO NOT USE
# 13187: CHD monitoring
# 13200: New patient health check
# 14448: Ambulatory blood pressure recording
# 14640: Lying blood pressure reading
# 14641: Standing blood pressure reading
# 14642: Standing blood pressure reading
# 14643: O/E - BP reading very high
# 15126: pre-treatment BP reading
# 18150: Coronary heart disease monitoring administration
# 18418: Sitting systolic blood pressure
# 19905: O/E - BP reading: no postural drop
# 20049: Blood pressure monitoring
# 22476: Blood pressure recorded by patient at home
# 22595: O/E - Diastolic BP reading; <- KEEP FOR DBP
# 23006: O/E - CVS examination NOS
# 23312: O/E - Systolic BP reading
# 25553: O/E - BP labile
# 27271: O/E - BP stable
# 27272: O/E - blood pressure reading NOS
# 27273: O/E - blood pressure decreased
# 27274: O/E - BP bordeline low
# 29261: Average home syslotic blood pressure
# 31305: Average home diastolic blood pressure <- KEEP FPR DBP
# 32208: O/E - Arterial pressure index abnormal
# 34618: Average day interval systolic blood pressure
# 34994: O/E - Arterial presure index normal
# 37242: Standing systolic blood pressure
# 37243: Standing diastolic blood pressure <- KEEP FOR DBP
# 38277: Average 24 hour systolic blood pressure
# 38278: Average night interval systolic blood pressure
# 41052: Ambulatory systolic blood pressure
# 41445: Average day interval diastolic blood pressure <- KEEP FOR DBP
# 42280: Sitting diastolic blood pressurel <- KEEP FOR DBP
# 43282: Lying systolic blood pressure
# 43547: Average night interval diastolic blood pressure <- KEEP FOR DBP
# 48008: Ambulatory diastolic blood pressure <- KEEP FOR DBP
# 65990: Lying diastolic blood pressure <- KEEP FOR DBP
# 94807: Self measured blood pressure reading
# 104203: Referral for ambulatory blood pressure monitoring

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
nrow(df_dbp)
# remove duplicates
df_dbp <- distinct(df_dbp)

# take averages for readings taken on the same day

df_dbp <- .dailymean(df_dbp, colname = "dbp")

# latest information
df_dbp <- .latest(df_dbp)

### add to df_base ###

df_base <- merge(df_base, df_dbp[, .(patid, dbp)], by = "patid", all.x = TRUE)
range(df_base$dbp, na.rm = TRUE)
sum(is.na(df_base$dbp)) / nrow(df_base)

######################################################################
### HDL cholesterol
######################################################################

# consider all tests first, as will be used later
df <- .clean_ctsvar(filename = "df_hdl.Rdata", colname = "hdl", pre = FALSE)

# entity types
filter(entity, enttype %in% unique(df$enttype))
# check medcodes
table(df$medcode)

# 44: Serum HDL cholesterol level
# 2379: seen in diabetic clinic; too few, DO NOT USE
# 13760: Serum fasting HDL cholesterol
# 13761: Serum random HDL cholesterol
# 13762: Serum plasma HDL cholesterol
# 14108: HDL:total cholesterol ratio, DO NOT USE
# 18147: Total cholesterol measurement
# 19635: Lipoprotein electroph. - HDL
# 26915: Plasma random HDL cholesterol
# 34548: Plasma HDL cholesterol level
# 93756: Non HDL cholesterol level, DO NOT USE
df <- df[!(medcode %in% c(2379, 11867, 93756))]

### units ###

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# usual units: mg/dl; mmol/L
df <- df[Specimen.Unit.Of.Measure %in% c("g/L", "mg/dL", "mmol/L", "mol/L")]

# convert into mmol/L
df[, hdl.mmol.L := NaN]

df[Specimen.Unit.Of.Measure == "g/L", hdl.mmol.L := hdl * 100 * 0.02586]
quantile(df$hdl.mmol.L[df$Specimen.Unit.Of.Measure == "g/L"]) 

df[Specimen.Unit.Of.Measure == "mg/dL", hdl.mmol.L := hdl * 0.02586]
quantile(df$hdl.mmol.L[df$Specimen.Unit.Of.Measure == "mg/dL"])

df[Specimen.Unit.Of.Measure == "mmol/L", hdl.mmol.L := hdl]
quantile(df$hdl.mmol.L[df$Specimen.Unit.Of.Measure == "mmol/L"])

df[Specimen.Unit.Of.Measure == "mol/L", hdl.mmol.L := hdl * 1000]
quantile(df$hdl.mmol.L[df$Specimen.Unit.Of.Measure == "mol/L"]) # too high; remove
df <- df[Specimen.Unit.Of.Measure != "mol/L"]

range(df$hdl.mmol.L)

# remove implausible values
df <- df[hdl.mmol.L >= 0.1 & hdl.mmol.L <= 30]

# mean per day
df <- .dailymean(df, colname = "hdl.mmol.L")

# save this copy to be used later
df_hdl <- df

# clean & add to the main dataset
df <- merge(df, df_studyentry, by = "patid", all.x = TRUE)
df <- df[eventdate <= studyentry]
df <- .latest(df)
df_base <- .add_ctsvar(df = df, colname = "hdl.mmol.L")
df_base <- rename(df_base, hdl = hdl.mmol.L)
df_base

######################################################################
### LDL cholesterol
######################################################################

# only tests taken before study entry are relevant
df <- .clean_ctsvar(filename = "df_ldl.Rdata", colname = "ldl", pre = TRUE)

# entity types
filter(entity, enttype %in% unique(df$enttype))
# check medcodes
table(df$medcode)

# 65: Serum LDL cholesterol level
# 13765: Serum fasting LDL cholesterol level
# 13766: Calculated LDL cholesterol level
# 19764: Plasma LDL Cholesterol level
# 29699: Plasma fasting LDL cholesterol level
# 33304: Plasma random LDL cholesterol level
# 46224: Serum random LDL cholesterol level
# 93756: Non HDL cholesterol level; DO NOT USE
# 108302: no medcode! DO NOT USE

df <- df[!(medcode %in% c(93756, 108302))]

### units ###

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# usual units: mg/dl; mmol/L
df <- df[Specimen.Unit.Of.Measure %in% c("g/L", "mg/dL", "mmol/L")]

# convert into mmol/L
df[, ldl.mmol.L := NaN]

df[Specimen.Unit.Of.Measure == "g/L", ldl.mmol.L := ldl * 100 * 0.02586]
quantile(df$ldl.mmol.L[df$Specimen.Unit.Of.Measure == "g/L"]) 

df[Specimen.Unit.Of.Measure == "mg/dL", ldl.mmol.L := ldl * 0.02586]
quantile(df$ldl.mmol.L[df$Specimen.Unit.Of.Measure == "mg/dL"])

df[Specimen.Unit.Of.Measure == "mmol/L", ldl.mmol.L := ldl]
quantile(df$ldl.mmol.L[df$Specimen.Unit.Of.Measure == "mmol/L"])

range(df$ldl.mmol.L)

# remove implausible values
df <- df[ldl.mmol.L >= 0.1 & ldl.mmol.L <= 30]

# clean & add to the main dataset
df_base <- .wrapper_ctsvar(df = df, colname = "ldl.mmol.L")
df_base <- rename(df_base, ldl = ldl.mmol.L)

######################################################################
### Total cholesterol
######################################################################

# all tests are needed, post-study entry tests will be used later
df <- .clean_ctsvar(filename = "df_serum_chol.Rdata", colname = "chol", pre = FALSE)

# entity types
filter(entity, enttype %in% unique(df$enttype))
# check medcodes
table(df$medcode)

# 0: no code; DO NOT USE
# 1: O/E - blood pressure reading; DO NOT USE
# 12: Serum cholesterol
# 44: Serum HDL cholesterol level; DO NOT USE
# 65: Serum LDL cholesterol level; DO NOT USE
# 241: Acute miocardial infarction; DO NOT USE
# 622: Serum cholesterol normal
# 2379: Seen in diabetic clinic; DO NOT USE
# 2493: Serum cholesterol raised
# 7041: Medication decreased; DO NOT USE
# 7913: Coronary heart disease risk; DO NOT USE
# 10940: Cholesterol screen
# 12821: Serum fasting total cholesterol
# 13733: Serum total cholesterol level
# 13760: Serum fasting HDL cholesterol level; DO NOT USE
# 13765: Serum fasting LDL cholesterol level; DO NOT USE
# 13770: Serum HDL:non-HDL cholesterol ratio; DO NOT USE
# 13771: Serum cholesterol/HDL ratio; DO NOT USE
# 13772: Total cholesterol:HDL ratio; DO NOT USE
# 14370: Serm HDL:non-HDL cholesterol ratio; DO NOT USE
# 14371: serum cholesterol/HDL ratio; DO NOT USE
# 14372: Total cholesterol : HDL ratio; DO NOT USE
# 18040: Plasma total cholesterol level
# 18147: Total cholesterol measurement
# 18443: pre-treatment serum cholesterol level
# 25180: Target cholesterol level; DO NOT USE
# 26902: Serum cholesterol NOS
# 29202: Serum cholesterol borderline 
# 34548: Plasma HDL cholesterol level; DO NOT USE
# 35720: Serum cholesterol very high
# 37206: Serum cholesterol studies; DO NOT USE
# 93756: Non HDL cholesterol level; DO NOT USE

df <- df[medcode %in% c(12, 622, 2493, 10940, 12821, 13733, 18040, 18147, 18443, 26902, 29202, 35720)]

### units ###

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# usual units: mg/dl; mmol/L
df <- df[Specimen.Unit.Of.Measure %in% c("g/dL", "mg/dL", "mg/L", "mmol/L", "mol/L", "umol/L")]
# NB: about 10% of data with "No Data Entered"; filtered out at the moment

# convert into mmol/L
df[, chol.mmol.L := NaN]

df[Specimen.Unit.Of.Measure == "g/dL", chol.mmol.L := chol * 1000 * 0.02586]
quantile(df$chol.mmol.L[df$Specimen.Unit.Of.Measure == "g/dL"]) 

df[Specimen.Unit.Of.Measure == "mg/dL", chol.mmol.L := chol * 0.02586]
quantile(df$chol.mmol.L[df$Specimen.Unit.Of.Measure == "mg/dL"])

df[Specimen.Unit.Of.Measure == "mg/L", chol.mmol.L := chol * 0.1 * 0.02586]
quantile(df$chol.mmol.L[df$Specimen.Unit.Of.Measure == "mg/L"])

df[Specimen.Unit.Of.Measure == "mmol/L", chol.mmol.L := chol]
quantile(df$chol.mmol.L[df$Specimen.Unit.Of.Measure == "mmol/L"])

df[Specimen.Unit.Of.Measure == "mol/L", chol.mmol.L := chol * 1000]
quantile(df$chol.mmol.L[df$Specimen.Unit.Of.Measure == "mol/L"]) # consistently too high; remove
df <- df[Specimen.Unit.Of.Measure != "mol/L"]

df[Specimen.Unit.Of.Measure == "umol/L", chol.mmol.L := chol * 0.001]
quantile(df$chol.mmol.L[df$Specimen.Unit.Of.Measure == "umol/L"]) # consistently too low; remove
df <- df[Specimen.Unit.Of.Measure != "umol/L"]

range(df$chol.mmol.L)

# remove implausible values
df <- df[chol.mmol.L >= 0.1 & chol.mmol.L <= 30]

# mean per day
df <- .dailymean(df, colname = "chol.mmol.L")

# save this copy to be used later
df_totchol <- df

# clean & add to the main dataset
df <- merge(df, df_studyentry, by = "patid", all.x = TRUE)
df <- df[eventdate <= studyentry]
df <- .latest(df)
df_base <- .add_ctsvar(df = df, colname = "chol.mmol.L")
df_base <- rename(df_base, chol = chol.mmol.L)
df_base
sum(is.na(df_base$chol.mmol.L)) / nrow(df_base)

######################################################################
### FH CKD
######################################################################

df <- .clean_df(filename = "df_fh_ckd.Rdata")
id <- unique(df$patid)
df_base[, fh_ckd := as.numeric(patid %in% id)]
table(df_base$fh_ckd, exclude = NULL)

######################################################################
### Kidney stones
######################################################################

#df <- .clean_df(filename = "df_kidney_stones.Rdata")
#id <- unique(df$patid)
#df_base[, fh_ckd := as.numeric(patid %in% id)]

######################################################################
### ACEIs
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_acei.Rdata")
df_base[, acei := as.numeric(patid %in% id_Tx)]
table(df_base$acei, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_acei.Rdata", days = 84)
df_base[, acei_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$acei, df_base$acei_84, exclude = NULL)

######################################################################
### ARBs
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_arb.Rdata")
df_base <- df_base[, arb := as.numeric(patid %in% id_Tx)]
table(df_base$arb, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_arb.Rdata", days = 84)
df_base[, arb_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$arb, df_base$arb_84, exclude = NULL)

######################################################################
### Diuretics
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_diur.Rdata")
df_base[, diur := as.numeric(patid %in% id_Tx)]
table(df_base$diur, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_diur.Rdata", days = 84)
df_base[, diur_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$diur, df_base$diur_84, exclude = NULL)

######################################################################
### NSAIDs
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_nsaid.Rdata")
df_base[, nsaid := as.numeric(patid %in% id_Tx)]
table(df_base$nsaid, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_nsaid.Rdata", days = 84)
df_base[, nsaid_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$nsaid, df_base$nsaid_84, exclude = NULL)

######################################################################
### Phosphate binders
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_pb.Rdata")
df_base[, pb := as.numeric(patid %in% id_Tx)]
table(df_base$pb, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_pb.Rdata", days = 84)
df_base[, pb_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$pb, df_base$pb_84, exclude = NULL)

######################################################################
### Beta blockers
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_beta.Rdata")
df_base[, beta := as.numeric(patid %in% id_Tx)]
table(df_base$beta, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_beta.Rdata", days = 84)
df_base[, beta_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$beta, df_base$beta_84, exclude = NULL)

######################################################################
### calcium channel blockers
######################################################################

# 28 days
id_Tx <-  .get_id_pre_Tx(filename = "df_ccb.Rdata")
df_base[, ccb := as.numeric(patid %in% id_Tx)]
table(df_base$ccb, exclude = NULL)

# 84 days
id_Tx_84 <-  .get_id_pre_Tx(filename = "df_ccb.Rdata", days = 84)
df_base[, ccb_84 := as.numeric(patid %in% id_Tx_84)]
table(df_base$ccb, df_base$ccb_84, exclude = NULL)

######################################################################
### Hypertension diagnosis 
######################################################################

### this is diagnosis only, without treatment

id_diag <- .get_id_pre_diag("df_hyp_diag.Rdata")
df_base[, hyp_diag := as.numeric(patid %in% id_diag)]
table(df_base$hyp_diag, exclude = NULL)
table(df_base$tr_hyp, df_base$hyp_diag, exclude = NULL)


######################################################################
### End of follow-up
######################################################################

### death date ###

# ONS: this is the most reliable source of information

df_1 <- df_death[, .(patid, deathdate)]

# df_pat: this is used only for identifying additional patients

df_pat <-  merge(df_pat, df_cens[, .(patid, end)])
df_pat
df_pat_death <- df_pat[!is.na(deathdate) & deathdate <= end]
id_1 <- df_1$patid
id_2 <- df_pat_death$patid
id <- setdiff(id_2, id_1) 
length(id) # 1,401 patients only
df_2 <- df_pat_death[patid %in% id]
df_2 <- df_2[, .(patid, deathdate)]
df_1 <- rbind(df_1, df_2)

# add to df_base
df_base <- merge(df_base, df_1, all.x = TRUE)
nrow(df_base)
range(df_base$deathdate, na.rm = TRUE)

### checks ###

### death date before study entry

alpha <- df_base[deathdate < studyentry]
alpha <- alpha[, .(patid, studyentry, deathdate)]
nrow(alpha) 
id <- alpha$patid
quantile(alpha$deathdate - alpha$studyentry)
# how many of these appear in the ONS?
beta <- df_death[patid %in% id]
nrow(beta) 
# assume date of study entry is incorrect, and therefore remove those patients
df_base <- df_base[!(patid %in% id)]
nrow(df_base)# NB: this will also remove patients whose CV event happened before study entry, as all of those came from ONS & were fatal

### end of follow-up ###

df_base[, fuenddate := pmin(deathdate, end, na.rm = TRUE)]
range(df_base$fuenddate) # everyone has an end of follow-up

## negative follow-up?

range(df_base$fuenddate - df_base$studyentry) # no, this has already been taken care of

table(df_base$diab_84)
df_base[is.na(diab_84), diab_84:=0]

######################################################################
### Total cholesterol : high density lipoprotein ratio
######################################################################

### Scenario 1: information on the ratio is available ###

#enttype_t <- subset(entity, description == "HDL/LDL ratio")$enttype # 338
# NB: this does NOT in fact correspond to what we want; but some entries do, so need to filter out carefully

# leave pre & post values
df <- .clean_ctsvar(filename = "df_chol.Rdata", colname = "ratio", pre = FALSE)

table(df$medcode, exclude = NULL)
# 44: Serum HDL cholesterol level; DO NOT USE
# 2379: Seen in diabetic clinic; DO NOT USE
# 14108: HDL: total cholesterol ratio; use but need to convert: x -> 1/x
# 14369: HDL : LDL ratio; DO NOT USE
# 14370: Serum HDL: non-HDL cholesterol ratio; use, but  need to convert: x -> 1/x + 1
# 14371: serum cholesterol/HDL ratio; use as is
# 14372: Total cholesterol : HDL ratio; use as is
# 19853: Serum LDL/HDL radio; DO NOT USE
# 35583: Serum cholesterol/LDL ratio: DO NOT USE
# 40935: Plasma cholsterol/HDL ratio; use as is
# 49695: Plasma LDL/HDL ratio; DO NOT USE
# 50393: Plasma cholesterol/LDL ratio; DO NOT USE
# 63314: Serum cholesterol/VLDL ratio; DO NOT USE
df <- df[medcode %in% c(14108, 14370, 14371, 14372, 40935)]

### units ###

table(df$Specimen.Unit.Of.Measure, exclude = NULL)
# usual units: 1/1, ratio
# Also leave "NO Data Entered" in this case, as this is the vast majority, and there really are no units
df <- df[Specimen.Unit.Of.Measure %in% c("1/1", "No Data Entered", "ratio")]

# convert into ratio; units are the same but numerators/denominators aren't

df[, chol_hdl := NaN]

df[medcode %in% c(14108), chol_hdl := 1 / ratio]
quantile(df$chol_hdl[df$medcode %in% c(14108)])
quantile(df$chol_hdl[df$medcode %in% c(14108)], probs = c(0.1, 0.5, 0.9, 0.95, 0.99)) # seems consistently <1, most likely recorded the other way round, but remove as not clear
# NB the ratio cannot be <1 by definition!
df <- df[medcode != 14108]

df[medcode %in% c(14370), chol_hdl := (1 / ratio) + 1]
quantile(df$chol_hdl[df$medcode %in% c(14370)])
quantile(df$chol_hdl[df$medcode %in% c(14370)], probs = c(0.1, 0.5, 0.9, 0.95, 0.95))
quantile(df$ratio[df$medcode %in% c(14370)], probs = c(0.1, 0.5, 0.9, 0.95, 0.95)) # whilst the values are >1 & hence plausible, it does seem that it's the total/hdl ratio that has been recorded; remove
df <- df[medcode != 14370]

df[medcode %in% c(14371, 14372, 40935), chol_hdl := ratio]
quantile(df$chol_hdl[df$medcode %in% c(14371, 14372, 40935)], probs = c(0.1, 0.5, 0.9, 0.95)) # just some few very abnormal values

# check
range(df$chol_hdl) 
subset(df, chol_hdl >= 50) # abnormally high
# also the values cannot be < 1 by definition
df <- df[chol_hdl >= 1 & chol_hdl <= 500] 

# clean & add to the output dataset

# take mean of tests taken on the same day
df <- .dailymean(df, colname = "chol_hdl")

df_ratio <- df

### Scenario 2: There are separate tests for Total & HDL cholesterol ###

# combine df_hdl & df_totchol
df <- merge(df_hdl, df_totchol, by = c("patid", "eventdate"))
df[, chol_hdl := chol.mmol.L / hdl.mmol.L]
df <- df[, .(patid, eventdate, chol_hdl)]

# check
quantile(df$chol_hdl)
df <- df[chol_hdl >= 1]

# combine with df_ratio
df <- rbind(df_ratio, df)
df <- distinct(df)

######################################################################
### FU dates, statins & incorporate Chol ratio
######################################################################

# add dates of study entry, cv and end of follow up
# end of follow up incorporates cv date where relevant
df <- merge(df, df_base[, .(patid, studyentry, fuenddate)], by = "patid")
df <- df[eventdate <= fuenddate]

# add date of statin prescription

# extract all precriptions
df_statin <- .clean_df(filename = "df_statins.Rdata")
# NB: at least two prescription would have been required at the extraction stage
# So just take the min date as the date of statin prescription
df_statin <- df_statin[order(patid, eventdate)]
df_statin <- df_statin[, head(.SD, 1), by = patid]
df_statin <- rename(df_statin, statindate = eventdate)

# add information to df
df <- merge(df, df_statin[, .(patid, statindate)], by = "patid", all.x = TRUE)
df[, censor := pmin(fuenddate, statindate, na.rm = TRUE)]
df <- df[eventdate <= censor]
df <- df[, .(patid, eventdate, chol_hdl)]

# extract value closest to  baseline
df <- .closest(df)

# add to df_base ###
df_base <- .add_ctsvar(df, colname = "chol_hdl")
sum(is.na(df_base$chol_hdl)) / nrow(df_base) 

######################################################################
### Further variables; added at a later stage
######################################################################

df_base_temp <- df_base

inputDir <- "K://SHARP_EE//Monitoring CKD CHF//Data//Derived CPRD data"

setwd(inputDir)

load("df_base_addvar.Rdata")

head(df_base)
new_drugs <- data.table(df_base)
new_drugs <- new_drugs[,.(patid, apc_base_84)]

df_base_temp <- merge(df_base_temp, new_drugs, by="patid", all.x=TRUE)
df_base_temp <- df_base_temp[!is.na(df_base_temp$apc_base_84),]
length(unique(df_base_temp$patid))
#df_base <- df_base_temp
setwd(sourceDir_R)

load("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\df_cancer_hf.Rdata")
df_cancer_hf <- df_base
df_cancer_hf <- df_cancer_hf[,.(patid, cancer)]

data <- merge(df_base_temp, df_cancer_hf, by="patid", all.x=TRUE)
data[,cancer:=ifelse(is.na(cancer),0,cancer)]

######################################################################
### Save
######################################################################

df_base <- data

save(df_base, file = paste0(outputFilename, ".Rdata"))