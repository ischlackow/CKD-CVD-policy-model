#############################################################
### Demographics
#############################################################

demographics <- function(df, vars, subgroupVar) {
  
  # vars is the list, each element of which has the format(list(varName = <colname in the dataframe>, cts = TRUE/FALSE, digits = integer)) 
  
  output <- NULL
  
  for (var in vars) {
    
    name <- var$varName
    digits <- var$digits
    
    if (var$cts) {
      
      ### if continuous variable, calculate mean and standard deviation
      
      # formula
      f <- as.formula(paste(subgroupVar, ".", sep = "~"))
      # numbers
      beta <- tapply(df[, name], df[, subgroupVar], function(x) paste(dp(mean(x, na.rm = TRUE), digits = digits), 
                                                                      " (", dp(sd(x, na.rm = TRUE), digits = digits), ")", sep = ""))
      # add overall
      beta <- c(var$varName, beta, all = paste(round(mean(df[, name], na.rm = TRUE), digits = digits), 
                                               " (", round(sd(df[, name], na.rm = TRUE), digits = digits), ")", sep = ""))
      
      # format & add to the output dataset
      names(beta)[length(names(beta))] <- "(all)"
      output <- if (is.null(output)) beta else rbind(output, beta) 
    } else {
      
      ### if discrete variable, calculate actual numbers and column percentages
      
      # formula
      f <- as.formula(paste(name, subgroupVar, sep = "~"))
      # numbers
      beta <- dcast(df, f, function(x) length(x), margins = subgroupVar)
      # add column percentages
      for (col in colnames(beta)[2 : ncol(beta)]){
        N <- sum(beta[, col])
        beta[, col] <- paste(beta[, col], " (", round(100 * beta[, col] / N, digits = digits), "%)", sep = "")
      }
      # format & add to the output dataset
      beta[, 1] <- as.character(beta[, 1])
      beta <- rbind(beta, c(colnames(beta)[1], rep(NA, ncol(beta) - 1)))
      beta <- beta[c(nrow(beta), 1 : (nrow(beta) - 1)), ]
      colnames(beta)[1] <- "var"
      
      output <- if (is.null(output)) beta else rbind(output, beta)
    }
  }
  return(data.table(output))
}

dp <- function(n, digits = 1) {
  
  # rounds, and displays, the number to the given number of decimal points
  formatC(round(n, digits = digits), digits = digits, format = "f")
  
}

.flag_inc <- function(output){
  output[, inc := 1]
  output[var %in% c("sex", "smok_short", "hyp_diag", "cvd", "diab", "cancer"), inc := 0]
  output[var == 0, inc := 0]
  output[var == "F", inc := 0]
  output[var %in% c("never", "before"), inc := 0]
  return(output)
}