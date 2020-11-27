######################################################################
### Functions
######################################################################

.clean_df <- function(filename){
  
  # load the file
  setwd(sourceDir_R_raw)
  df <- data.table(get(load(filename)))
  
  # leave only entities corresponding to numeric test result
  enttype_t <- entity$enttype[entity$enttype %in% unique(df$enttype) & entity$data1 == "Operator"]
  df <- df[enttype %in% enttype_t]
  df[, enttype := NULL]
  
  ### clean the data
  
  df <- df[, .(patid, eventdate, data1, data2, data3)]
  # always remove NA
  df <- df[!is.na(data2)]
  #if (remove_zero){
  #  # remove 0 values as these are either not interpretable or a non-feasible value anyway
  #  df <- df[data2 > 0]
  #}
  df <- distinct(df)
  
  # operators
  df <- rename(df, Code = data1)
  df <- merge(df, OPR, by = "Code")
  df[, Code := NULL]
  df[, Operator := drop.levels(Operator)]
  # remove "Data Not Entered" operator unless the value is 0
  df <- df[(data2 != 0 & Operator != "Data Not Entered") | (data2 == 0)]
  
  # units
  df <- rename(df, Code = data3)
  df <- merge(df, SUM, by = "Code")
  df[, Code := NULL]
  df[, Specimen.Unit.Of.Measure := drop.levels(Specimen.Unit.Of.Measure)]
  
  return(df)
  
}

.filter_lb <- function(df, colname, lb){
  df <- df[(Operator %in% c("=", "~")) | 
             (Operator %in% c("<", "<=") & get(colname) < lb) |
             (Operator %in% c("<") & get(colname) == lb) |
             (Operator == "Data Not Entered" & get(colname) == 0)]
  return(df)
}

.filter_ub <- function(df, colname, ub){
  df <- df[(Operator %in% c("=", "~")) | 
             (Operator %in% c(">", ">=") & get(colname) > ub) |
             (Operator %in% c(">") & get(colname) == ub) |
             (Operator == "Data Not Entered" & get(colname) == 0)]
  return(df)
}


.filter_lb_up <- function(df, colname, ub, lb){
  df <- df[(Operator %in% c("=", "~")) | 
               (Operator %in% c("<", "<=") & get(colname) < lb) |
               (Operator %in% c("<") & get(colname) == lb) | 
               (Operator %in% c(">", ">=") & get(colname) > ub) |
               (Operator %in% c(">") & get(colname) == ub) |
               (Operator == "Data Not Entered" & get(colname) == 0)]
  return(df)
}



.clean_df_post <- function(df, colname, mean_per_day = TRUE){
  
  # remove redundant columns
  df[, Specimen.Unit.Of.Measure := NULL]
  df[, data2 := NULL]
  
  # mean per day
  if (mean_per_day)
    df <- df[, lapply(.SD, mean), by = .(patid, eventdate), .SDcols = colname]
  
  return(df)
}

.assign_stage_cts <- function(df, lb, ub) {
  
  df[, stage := NaN]
  df[val < lb, stage := 1]
  df[val >= lb & val <= ub, stage := 2]
  df[val > ub, stage := 3]
  
  return(df)
}