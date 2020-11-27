### combine 4 simulation datasets into 1 ###########################################

.combine_4_into_1 <- function(file_prefix, file_suffix = NULL, dir) {
  
  setwd(dir)
  
  df <- NULL
  for (sex in c("F", "M"))
    for (secondary in c(0, 1)) {
      temp <- get(load(paste0(file_prefix, "_", sex, "_", secondary, file_suffix, ".Rdata")))
      temp <- mutate(temp, secondary = secondary, sex = sex)
      df <- rbind(df, temp)
    }
  
  return(df)
}

### cumulative incidence ###########################################

.cumprob_Act <- function(dfAct, Y, subgroups, endpts) {
  
  output <- NULL
  
  for (endpt in endpts) {
    
    # read off endpt information
    lab <- endpt$lab
    endptVar <- endpt$endptAct
    fuVar <- endpt$fuvarAct
    
    # events per year
    for (subgroup in subgroups) {
      
      dfTemp <- if (subgroup == "all") dfAct else filter(dfAct, ckd == subgroup)
      if (nrow(dfTemp) > 0) {
        # fit Kaplan-Meier estimate
        f <- as.formula(paste("Surv(", fuVar, ", ", endptVar, ") ~ 1", sep = ""))
        km <- survfit(f, data = dfTemp)
        risk <- data.frame(time = km$time, cumprob = km$surv, lower = km$lower, upper = km$upper)
        risk <- mutate(risk, year = ceiling(time))
        alpha <- risk %>%
          group_by(year) %>%
          summarise(est = 1 - min(cumprob), l = 1 - min(upper), u = 1 - min(lower))  %>% # NB upper & lower need to be swapped as negation
          ungroup() %>% 
          filter(year <= Y & year > 0) %>%
          arrange(year)
        if (nrow(alpha) > 0) 
          alpha <- alpha %>%
          mutate(subgroup = subgroup, lab = lab)
        output <- rbind(output, alpha)
      }
    }
  }
  return(output)
}

.cumprob_Pred <- function(dfPred, subgroups, endpts, Y) {
  
  output <- NULL
  
  for(endpt in endpts) {
    
    # read off endpt information
    lab <- endpt$lab
    endptVar <- endpt$endptPred 
    fuVar <- endpt$fuvarPred
    #t <- 1
    # events per year
    for (t in 1 : Y) {
      dfTemp <- filter(dfPred, cycle == t)
      for (subgroup in subgroups) {
        dfTemp2 <- if (subgroup == "all") dfTemp else filter(dfTemp, ckd == subgroup)
        if (nrow(dfTemp2) > 0) {
          N_events <- sum(dfTemp2[, endptVar])
          N_risk <- sum(dfTemp2[, fuVar])
          output <- rbind(output, data.frame(lab = lab, year = t, subgroup = subgroup, prod = 1-N_events/N_risk))
        }
      }
    }
  }
  # KM product
  output <- output %>%
    group_by(lab, subgroup) %>%
    mutate(km = 1-cumprod(prod))
  return(output)
}

### graphical ###########################################

.get_p_intval <- function(df, year, p_breaks, x_axis_label, y_axis_label) {
  
  p_labels <- paste0(p_breaks * 100, "%")
  
  p <- ggplot(df, aes(x = year))
  # actual
  p <- p + geom_errorbar(aes(x = year - 0.15, ymin = l, ymax = u), width = 0.35, size = 2)
  p <- p + geom_point(aes(y = est, x = year - 0.15), size = 4)
  # predicted
  p <- p + geom_point(aes(x = year + 0.15, y = km), size = 4, colour = "#665D5D")
  # format
  p <- p + facet_grid(lab ~ subgroup2, labeller = labeller(type = label_parsed))
  p <- p + coord_cartesian(xlim = c(0.7, Y + 0.6))
  p <- p + scale_x_continuous(name = x_axis_label, breaks = c(1 : Y))
  p <- p + scale_y_continuous(name = y_axis_label, breaks = p_breaks, labels = p_labels)
  p <- p + theme_bw() + theme(legend.key = element_blank(), 
                              axis.text = element_text(size = 23), 
                              axis.title = element_text(size = 25),
                              strip.text.x = element_text(size = 25),
                              strip.text.y = element_text(size = 25),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor = element_blank())
  
  return(p)
}

### collate PSA results ###########################################

### applications ###########################################

## data preparation for the application
.clean_df_app <- function(dir, file_prefix, file_suffix = NULL){
  
  # load & combine the files
  df <- .combine_4_into_1(file_prefix = file_prefix, file_suffix = file_suffix, dir = dir)
  
  ### add baseline info
  setwd(data_for_extrapolation_dir)
  df_base <- get(load("df_sim_id.Rdata"))
  df_base <- select(df_base, patid_temp, age_cat, ckd)
  df <- merge(df, df_base)
  
  df <- df %>% 
    select(-sex) %>%
    select(-patid_temp) %>%
    group_by_at(vars) %>%
    summarise_each(mean)
  
  return(df)
}

.clean_df_psa <- function(dir, file_prefix, n) {
  
  df <- .clean_df_app(dir = dir, file_prefix = file_prefix, file_suffix = paste0("_", n))
  df <- df %>% 
    mutate(sim = n) %>%
    select(sim, everything())
  
  return(df)
}

ci_psa <- function(v, p) {
  v <- sort(v)
  n <- length(v)
  ll <- ceiling(n * p / 2)
  ul <- floor(n * (1 - p / 2))
  return(list(l = v[ll], u = v[ul]))
}

.collate_df_psa <- function(df, vars, p = 0.05) {
  
  # lower limit
  df_l <- df %>%
    ungroup() %>%
    group_by_at(vars) %>%
    summarise(ly = ci_psa(ly, p)[["l"]], qaly = ci_psa(qaly, p)[["l"]]) %>%
    mutate(limit = "lower")
  
  df_u <- df %>%
    ungroup() %>%
    group_by_at(vars) %>%
    summarise(ly = ci_psa(ly, p)[["u"]], qaly = ci_psa(qaly, p)[["u"]]) %>%
    mutate(limit = "upper")
  
  df <- rbind(df_l, df_u)
  
  return(df)
  
}

.clean_df_undertreat <- function(dir_0, file_prefix_0, dir_1 = dir_0, file_prefix_1, file_suffix = NULL, p_undertreat){
  
  df_0 <- .clean_df_app(dir = dir_0, file_prefix = file_prefix_0, file_suffix = file_suffix)
  df_1 <- .clean_df_app(dir = dir_1, file_prefix = file_prefix_1, file_suffix = file_suffix)
  df <- merge(df_0, df_1, by = vars, suffix = c("_0", "_1"))
  
  # add undertreatment value
  df <- df %>%
    mutate(ly_p = ly_0 * (1 - p_undertreat[secondary + 1]) + ly_1 * p_undertreat[secondary + 1],
           qaly_p = qaly_0 * (1 - p_undertreat[secondary + 1]) + qaly_1 * p_undertreat[secondary + 1]) %>%
    select(secondary, ckd, age_cat, everything()) %>%
    mutate(ly_inc_0_p = ly_p - ly_0, qaly_inc_0_p = qaly_p - qaly_0,
           ly_inc_p_1 = ly_1 - ly_p, qaly_inc_p_1 = qaly_1 - qaly_p) %>%
    select(-(ly_0:qaly_p))
  
  return(df)
  
}

.clean_df_undertreat_psa <- function(dir_0, file_prefix_0, dir_1, file_prefix_1, p_undertreat, n) {
  
  df <- .clean_df_undertreat(dir_0 = dir_0, file_prefix_0 = file_prefix_0, 
                             dir_1 = dir_1, file_prefix_1 = file_prefix_1, 
                             file_suffix =  paste0("_", n), 
                             p_undertreat = p_undertreat)
  df <- df %>% 
    mutate(sim = n) %>%
    select(sim, everything())
  
  return(df)
}

.collate_df_undertreat_psa <- function(df, vars, p = 0.05) {
  
  # lower limit
  df_l <- df %>%
    ungroup() %>%
    group_by_at(vars) %>%
    summarise(ly_inc_0_p = ci_psa(ly_inc_0_p, p)[["l"]], qaly_inc_0_p = ci_psa(qaly_inc_0_p, p)[["l"]],
              ly_inc_p_1 = ci_psa(ly_inc_p_1, p)[["l"]], qaly_inc_p_1 = ci_psa(qaly_inc_p_1, p)[["l"]]) %>%
    mutate(limit = "lower")
  
  df_u <- df %>%
    ungroup() %>%
    group_by_at(vars) %>%
    summarise(ly_inc_0_p = ci_psa(ly_inc_0_p, p)[["u"]], qaly_inc_0_p = ci_psa(qaly_inc_0_p, p)[["u"]],
              ly_inc_p_1 = ci_psa(ly_inc_p_1, p)[["u"]], qaly_inc_p_1 = ci_psa(qaly_inc_p_1, p)[["u"]]) %>%
    mutate(limit = "upper")
  
  df <- rbind(df_l, df_u)
  
  return(df)
  
}

.get_p_ly <- function(df, title, xvar, yvar, adj, ylims, ylims2, y_title, 
                      top_graph = FALSE, left_graph = FALSE, bottom_graph = FALSE, right_graph = FALSE){
  
  df <- mutate(df,label1 = round(value, 2), ylabel1 = value + adj)
  
  p1 <- ggplot(df, aes_string(x = xvar))
  p1 <- p1 + geom_bar(stat = "identity", aes(y = value), fill = "#B0AAAA", colour = "#555252")
  p1 <- p1 + facet_grid(. ~ ckd)
  p1 <- p1 + scale_x_discrete(name = "Age at cohort entry, years") 
  p1 <- p1 + theme_bw() + theme(
    axis.text.x = element_text(size = 17, angle = 60, hjust=0.95, colour = "black"),
    axis.text.y = element_text(size = 17),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 21),
    #strip.background = element_blank(),
    axis.title.x = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.7),
    plot.title = element_text(size = 27, hjust = 0.5))
  p1 <- p1 + scale_y_continuous(name = y_title, limits = ylims2, breaks = ylims, expand = c(0, 0))
  # add numbers
  p1 <- p1 + geom_text(df, mapping = aes(label = sprintf('%.1f',df$label1), x = data.frame(df)[, xvar] , y = ylabel1), 
                       size = 7, fontface = "bold")
  p1 <- p1 + ggtitle(title)
  if (top_graph)
    p1 <- p1 + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  if (bottom_graph)
    p1 <- p1 + theme(
      strip.text = element_blank(),
      plot.title = element_blank()
    )
  if (right_graph)
    p1 <- p1 + theme(
      axis.title.y = element_blank()
    )
  return(p1)
}

.get_p_ly_undertreat <- function(df, title, xvar, yvar1, yvar2, adj, ylims, ylims2, y_title, 
                                 top_graph = FALSE, left_graph = FALSE, bottom_graph = FALSE, right_graph = FALSE){
  
  df$label1 <- round(df[, yvar1], 2)
  df$label2 <- round(df[, yvar1] + df[, yvar2], 2)
  df$ylabel1<-  df[, yvar1] + adj
  df$tot <- df[, yvar1] + df[, yvar2]
  df$ylabel2 <- df$tot + adj
  
  p1 <- ggplot(df, aes_string(x = xvar))
  p1 <- p1 + geom_bar(stat = "identity", aes_string(y = "tot"), fill = "#B0AAAA", colour = "#555252")
  p1 <- p1 + geom_bar(stat = "identity", aes_string(y = yvar1), fill = "#555252", colour = "#555252")
  p1 <- p1 + facet_grid(. ~ ckd)
  p1 <- p1 + scale_x_discrete(name="Age at cohort entry, years") 
  p1 <- p1 + scale_y_continuous(name = y_title, limits = ylims2, breaks = ylims, expand = c(0, 0))
  ## add numbers
  #p1 <- p1 + geom_text(df, mapping = aes(label = sprintf('%.1f',df$label1), x = data.frame(df)[,xvar] , y = ylabel1), 
  #                     size = 5, colour = "white", fontface = "bold")
  #p1 <- p1 + geom_text(df, mapping = aes(label = sprintf('%.1f',df$label2), x = data.frame(df)[,xvar] , y = ylabel2), 
  #                     size = 5, fontface = "bold")
  p1 <- p1 + ggtitle(title)
  p1 <- p1 + theme_bw() + theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust=0.95),
    axis.text.y = element_text(size = 13),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 15),
    #strip.background = element_blank(),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.7),
    plot.title = element_text(size = 23, hjust = 0.5))
  if (top_graph)
    p1 <- p1 + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  if (bottom_graph)
    p1 <- p1 + theme(
      strip.text = element_blank(),
      plot.title = element_blank()
    )
  if (right_graph)
    p1 <- p1 + theme(
      axis.title.y = element_blank()
    )
  return(p1)
}
