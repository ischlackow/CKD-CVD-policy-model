### data extraction from stata files ###

# re-write
.getcodes <- function(sourceDir_codelist, filename, colname){
  
  setwd(sourceDir_codelist)
  codes_t <- fread(filename, header = T)
  
  return(codes_t$colname[included == 1])
  
}

# re-write
.getid_medcodes <- function(sourceDir_codelist, filename, df_clinical_t, df_ref_t){
  
  medcodes_t <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = filename, colname = "medcode")
  
  id_clinical <- unique(filter(df_clinical_t, medcode %in% medcodes_t)$patid)
  id_ref <- unique(filter(df_ref_t, medcode %in% medcodes_t)$patid)
  
  return(union(id_clinical, id_ref))
}

# re-write
.getid_Tx <- function(sourceDir_codelist, filename, df_therapy_t, days = 28){
  
  prodcodes_t <- .getcodes(sourceDir_codelist = sourceDir_codelist, filename = filename, colname = "prodcode")
  
  # extract all prescriptions
  df <- filter(df_therapy_t, prodcode %in% prodcodes_t)
  df <- select(df, patid, eventdate, studyentry)
  df <- distinct(df)
  
  # require at least two prescriptions
  df <- group_by(df, patid)
  alpha <- summarise(df, n = length(eventdate))
  beta <- filter(alpha, n >= 2)
  
  # of these, require the latest prescription to be at most <days> days away from studyentry
  id <- beta$patid
  df <- filter(df, patid %in% patid)
  df <- mutate(df, diff = as.numeric(studyentry - eventdate))
  df <- filter(df, diff <= days)
  df <- ungroup(df)
  
  return(unique(df$patid))
}

# write & re-write
.getid_test <- function(sourceDir_codelist, filename){
  
  
  return(id)
}

### data manipulation ###

.clean_df <- function(sourceDir = sourceDir_R_raw, filename, eventdate_colname = "eventdate", cens = FALSE){
  
  setwd(sourceDir)
  
  df <- data.table(get(load(filename)))
  df <- df[patid %in% id_ckd]
  
  if (cens) {
    df <- df[df_cens, nomatch = 0, on = "patid"]
    df <- df[get(eventdate_colname) >= start & get(eventdate_colname) <= end]
  }
  return(df)
}

.clean_df_all <- function(sourceDir = sourceDir_R_raw, filename, eventdate_colname = "eventdate", cens = FALSE){
  
  setwd(sourceDir)
  
  df <- data.table(get(load(filename)))
  #df <- df[patid %in% id_ckd]
  
  if (cens) {
    df <- df[df_cens, nomatch = 0, on = "patid"]
    df <- df[get(eventdate_colname) >= start & get(eventdate_colname) <= end]
  }
  return(df)
}

.split_df <- function(df, eventdate_colname = "eventdate") {
  
  df <- df[df_studyentry, on = "patid"]
  
  dfpre <- df[get(eventdate_colname) <= studyentry]
  dfpost <- df[get(eventdate_colname) > studyentry]
  
  return(list(dfpre = dfpre, dfpost = dfpost))
}

.get_id_pre_diag <- function(filename, cens = FALSE) {
  
  df <- .clean_df(filename = filename, cens = cens)
  df <- .split_df(df)$dfpre
  
  return(unique(df$patid))
}

.get_id_pre_Tx <- function(filename, cens = FALSE, days = 28) {
  
  setwd(sourceDir_R_raw)
  
  # extract all precriptions
  df <- .clean_df(filename = filename, cens = cens)
  df <- .split_df(df)$dfpre
  
  # require at least two prescriptions
  alpha <- df[, .N , by = patid]
  alpha <- alpha[N >= 2]
  
  # of these, require the latest prescription to be at most <days> days away from studyentry
  # default is 28, but also generate for 84 days
  id <- alpha$patid
  df <- df[patid %in% id]
  df <- df[, diff := as.numeric(studyentry - eventdate)]
  df <- df[diff <= days]
  
  return(unique(df$patid))
}

.get_id_Tx_pred1 <- function(filename, d1, cens = FALSE, days = 28) {
  
  setwd(sourceDir_R_raw)
  
  # extract all precriptions by date d1
  df <- .clean_df(filename = filename, cens = cens)
  df <- df[eventdate <= d1]
  
  # require at least two prescriptions
  alpha <- df[, .N , by = patid]
  alpha <- alpha[N >= 2]
  
  # of these, require the latest prescription to be at most <days> days away from studyentry
  # default is 28, but also generate for 84 days
  id <- alpha$patid
  df <- df[patid %in% id]
  df <- df[, diff := as.numeric(d1 - eventdate)]
  df <- df[diff <= days]
  
  return(unique(df$patid))
}

.get_id_Tx_d1d2 <- function(filename, d1, d2, cens = FALSE, days = 28) {
  
  setwd(sourceDir_R_raw)
  
  # extract all precriptions by date d1
  df <- .clean_df(filename = filename, cens = cens)
  df <- df[eventdate >= d1 & eventdate <= d2]
  
  # require at least two prescriptions
  alpha <- df[, .N , by = patid]
  alpha <- alpha[N >= 2]
  
  # no further restrictions
  
  return(unique(alpha$patid))
}

.get_id_preFU_Tx <- function(filename, cens = FALSE, days = 28) {
  
  setwd(sourceDir_R_raw)
  
  # extract all precriptions
  df <- .clean_df(filename = filename, cens = cens)
  df <- merge(df, df_FU, by = "patid") # do not split; doesn't matter whether the prescription was at studyentry or prior
  
  # require at least two prescriptions
  alpha <- df[, .N , by = patid]
  alpha <- alpha[N >= 2]
  
  # of these, require the latest prescription to be at most <days> days away from fuenddate
  # default is 28, but also generate for 84 days
  id <- alpha$patid
  df <- df[patid %in% id]
  df <- df[, diff := as.numeric(fuenddate - eventdate)]
  df <- df[diff <= days]
  df <- df
  
  return(unique(df$patid))
}

### clean a continuous variable ###

.add_units <- function(df) {
  df <- rename(df, Code = data3)
  df <- df[SUM, nomatch = 0, on = "Code"]
  df <- df[, Code:= NULL]
  df[, Specimen.Unit.Of.Measure := drop.levels(Specimen.Unit.Of.Measure)]
  
  return(df)
}

.add_operator <- function(df, oper_to_leave = "=") {
  
  df <- rename(df, Code = data1)
  df <- df[OPR, nomatch = 0, on = "Code"]
  df <- df[, Code:= NULL]
  df[, Operator := drop.levels(Operator)]
  
  if (oper_to_leave != FALSE)
    df <- df[Operator %in% oper_to_leave]
  
  return(df)
}

.clean_ctsvar <- function(filename, colname, oper_to_leave = "=", pre = TRUE){
  
  # extract sensible data
  df <- .clean_df(filename = filename)
  if (pre)
    df <- .split_df(df)$dfpre
  df <- df[!is.na(data2)]
  df <- distinct(df)
  
  # rename data2 column
  df[, eval(colname) := data2]
  df[, data2 := NULL]
  
  # add units & filter by operator
  df <- .add_units(df)
  df <- .add_operator(df, oper_to_leave)
  
  return(df)
}

.dailymean <- function(df, colname){
  
  df <- df[, lapply(.SD, mean), by = .(patid, eventdate), .SDcols = colname]
  
  return(df)
  
}

.latest <- function(df){
  
  df <- df[order(patid, eventdate)]
  df <- df[, tail(.SD, 1), by = patid]
  
  return(df)
}

.closest <- function(df){
  
  df <- df_base[, .(patid, studyentry)][df, nomatch = 0, on = "patid"]
  #df <- df[, .(patid, studyentry, eventdate, category)]
  df[, diff := as.numeric(abs(eventdate - studyentry))]
  
  df <- df[order(patid, diff)]
  df <- df[, head(.SD, 1), by = patid]
  
  return(df)
  
}

.add_ctsvar <- function(df, colname){
  
  df <- df[, c("patid", colname), with = FALSE]
  df_base <- merge(df_base, df, by = "patid", all.x = TRUE)
  
  return(df_base)
}

.wrapper_ctsvar <- function(df, colname){
  
  # take mean of tests taken on the same day
  df <- .dailymean(df, colname = colname)
  
  # pick the latest test
  df <- .latest(df)
  
  # add to the main dataset
  df_base <- .add_ctsvar(df, colname = colname)
  
  return(df_base)
}