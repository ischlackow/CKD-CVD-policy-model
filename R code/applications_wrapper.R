
### life expectancy assuming no treatment ----------------------------

rm(list = ls())

library(here)
library(tidyverse)

### directories & functions

source(here("src", "config", "define_dir.R"))
source(here("src", "model_paper", "model_paper_functions.R"))

### parameters

vars <- c("secondary", "ckd", "age_cat")

file_prefix <- "df_output_no_treat"

### extract data: point estimates

df_est <- .clean_df_app(dir = datasets_output_dir, file_prefix = file_prefix)

### extract data: psa

library(parallel)
library(doParallel)

n_sim <- 1000

input_dir <- file.path(datasets_output_dir, "psa", "app", "no_treat")

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, library(tidyverse))
retval <- foreach (n = 1:n_sim) %dopar%
  return(.clean_df_psa(dir = input_dir, file_prefix = file_prefix, n = n))
stopCluster(cl)
df_psa <- as.data.frame(do.call("rbind", retval))

### table

# point estimates
library(reshape2)

alpha <- melt(df_est, id.vars = vars)
beta <- dcast(alpha, variable + ckd ~ secondary + age_cat, value.var = "value")
beta

# uncertainty

p <- 0.05
df_psa_sum <- .collate_df_psa(df = df_psa, vars = vars, p = p) 
alpha_psa <- melt(df_psa_sum, id.vars = c(vars, "limit"))
beta_psa <- dcast(alpha_psa, variable + ckd ~ secondary + age_cat + limit, value.var = "value")
beta_psa <- select(beta_psa, variable, ckd, ends_with("_lower"), ends_with("_upper"))

# combine & save

output <- merge(beta, beta_psa, by = c("variable", "ckd"), sort = FALSE)
#output

setwd(model_paper_dir)
write.csv(output, file = "ly_qaly_no_treat.csv", row.names = FALSE)

### graph

library(ggplot2)
library(gridExtra)
library(grid)

ckd_labels <- c(
  paste0("eGFR 60", "\U2013", "89", "\n mL/min/1.73", "\U00B2", "\n\n(G2)"),
  paste0("eGFR 45", "\U2013", "59", "\n mL/min/1.73", "\U00B2", "\n\n(G3a)"),
  paste0("eGFR 30", "\U2013", "44", "\n mL/min/1.73", "\U00B2", "\n\n(G3b)"),
  paste0("eGFR 15", "\U2013", "29", "\n mL/min/1.73", "\U00B2", "\n\n(G4)"),
  paste0("eGFR <15", "\nmL/min/1.73", "\U00B2", "\nnot on RRT\n(G5)")
)

age_labels <- c("<60", "60-69", "70-79", "80+")

beta <- alpha %>%
  ungroup() %>%
  mutate(ckd = as.numeric(ckd), age = as.numeric(age_cat)) %>%
  mutate(ckd = factor(ckd, levels = 1:length(ckd_labels), labels = ckd_labels),
         age = factor(age, levels = 1:length(age_labels), labels = age_labels)) %>%
  select(-age_cat)


title_no_cv <- "Without cardiovascular disease\n"
title_cv <- "With cardiovascular disease\n"
xvar <- "age"

# life-years
yvar <- "ly"
adj <- 1.5
ylims <- seq(10, 40, by = 10)
ylims2 <- c(0, 42)
y_title <- "Survival, years"
p_ly_prim <- .get_p_ly(df = filter(beta, secondary == 0 & variable == yvar), 
                       title = title_no_cv, xvar = xvar, adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title,
                       top_graph = TRUE, left_graph = TRUE)
p_ly_sec <- .get_p_ly(df = filter(beta, secondary == 1 & variable == yvar), 
                      title = title_cv, xvar = xvar, adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title,
                      top_graph = TRUE, right_graph = TRUE)

# qalys
yvar <- "qaly"
adj <- 1.5
ylims <- seq(10, 30, by = 10)
ylims2 <- c(0, 35)
y_title <- "Quality adjusted survival, QALYs"
p_qaly_prim <- .get_p_ly(df = filter(beta, secondary == 0 & variable == yvar), 
                       title = title_no_cv, xvar = xvar, adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title,
                       bottom_graph = TRUE, left_graph = TRUE)
p_qaly_sec <- .get_p_ly(df = filter(beta, secondary == 1 & variable == yvar), 
                      title = title_cv, xvar = xvar, adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title,
                      bottom_graph = TRUE, right_graph = TRUE)

# combine
setwd(model_paper_dir)
png("ly_qaly_no_treat.png", height = 150+130, width = 220+200, units = "mm", res = 1000) 
grid.newpage()
pushViewport(viewport(
  height = unit(150+130, "mm"), width = unit(220+200, "mm"), 
  layout = grid.layout(2, 2, heights = unit(c(150, 130), "mm"), widths = unit(c(220, 200), "mm"))))
print(p_ly_prim, 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_ly_sec, 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_qaly_prim, 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_qaly_sec, 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()


### impact of undertreatment --------------------------

rm(list = ls())

library(here)
library(tidyverse)

### directories & functions

source(here("src", "config", "define_dir.R"))
source(here("src", "model_paper", "model_paper_functions.R"))

### parameters

vars <- c("secondary", "ckd", "age_cat")
p_undertreat <- c(prim = 0.37, sec = 0.71)

file_prefix_0 <- "df_output_no_treat"
file_prefix_1 <- "df_output_full_treat"

### extract data: point estimates
df_est <- .clean_df_undertreat(dir_0 = datasets_output_dir, 
                               file_prefix_0 = file_prefix_0, file_prefix_1 = file_prefix_1, file_suffix = NULL, 
                               p_undertreat = p_undertreat)

### extract data: psa

library(parallel)
library(doParallel)

n_sim <- 1000

input_dir_0 <- file.path(datasets_output_dir, "psa", "app", "no_treat")
input_dir_1 <- file.path(datasets_output_dir, "psa", "app", "full_treat")

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, library(tidyverse))
retval <- foreach (n = 1:n_sim) %dopar%
  return(.clean_df_undertreat_psa(dir_0 = input_dir_0, file_prefix_0 = file_prefix_0, 
                                  dir_1 = input_dir_1, file_prefix_1 = file_prefix_1, 
                                  p_undertreat = p_undertreat, n = n)) 
stopCluster(cl)
df_psa <- as.data.frame(do.call("rbind", retval))

### table

# point estimates
library(reshape2)

alpha <- melt(df_est, id.vars = vars)
beta <- dcast(alpha, variable + ckd ~ secondary + age_cat, value.var = "value")
beta[, c(1, 2, 4)]
beta[, c(1, 2, 9)]

# uncertainty

p <- 0.05
df_psa_sum <- .collate_df_undertreat_psa(df = df_psa, vars = vars, p = p) 
alpha_psa <- melt(df_psa_sum, id.vars = c(vars, "limit"))
beta_psa <- dcast(alpha_psa, variable + ckd ~ secondary + age_cat + limit, value.var = "value")
beta_psa <- select(beta_psa, variable, ckd, ends_with("_lower"), ends_with("_upper"))

# combine & save

output <- merge(beta, beta_psa, by = c("variable", "ckd"), sort = FALSE)
#output

setwd(model_paper_dir)
write.csv(output, file = "ly_qaly_inc.csv", row.names = FALSE)

### graph
  
library(ggplot2)
library(gridExtra)
library(grid)

ckd_labels <- c(
  paste0("eGFR 60", "\U2013", "89", "\n mL/min/1.73", "\U00B2", "\n\n(G2)"),
  paste0("eGFR 45", "\U2013", "59", "\n mL/min/1.73", "\U00B2", "\n\n(G3a)"),
  paste0("eGFR 30", "\U2013", "44", "\n mL/min/1.73", "\U00B2", "\n\n(G3b)"),
  paste0("eGFR 15", "\U2013", "29", "\n mL/min/1.73", "\U00B2", "\n\n(G4)"),
  paste0("eGFR <15", "\nmL/min/1.73", "\U00B2", "\nnot on RRT\n(G5)")
)

age_labels <- c("<60", "60-69", "70-79", "80+")

beta <- alpha %>%
  ungroup() %>%
  mutate(ckd = as.numeric(ckd), age = as.numeric(age_cat)) %>%
  mutate(ckd = factor(ckd, levels = 1:length(ckd_labels), labels = ckd_labels),
         age = factor(age, levels = 1:length(age_labels), labels = age_labels)) %>%
  select(-age_cat)

gamma <- dcast(beta, secondary + ckd + age ~ variable, value.var = "value")

title_no_cv <- "Without cardiovascular disease\n"
title_cv <- "With cardiovascular disease\n"
xvar <- "age"

# life-years
yvar1 <- "ly_inc_0_p"
yvar2 <- "ly_inc_p_1"
adj <- -0.05
ylims <- seq(0.5, 2.0, by = 0.5)
ylims2 <- c(0, 2.2)
y_title <- "Additional survival, \nyears"
p_ly_prim <- .get_p_ly_undertreat(df = filter(gamma, secondary == 0), 
                                  title = title_no_cv, 
                                  xvar = xvar, yvar1 = yvar1, yvar2 = yvar2, 
                                  adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title, 
                                  top_graph = TRUE, left_graph = TRUE)
p_ly_sec <- .get_p_ly_undertreat(df = filter(gamma, secondary == 1), 
                                  title = title_cv, 
                                  xvar = xvar, yvar1 = yvar1, yvar2 = yvar2, 
                                  adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title, 
                                  top_graph = TRUE, right_graph = TRUE)

# qalys
yvar1 <- "qaly_inc_0_p"
yvar2 <- "qaly_inc_p_1"
adj <- -0.05
ylims <- seq(0.5, 2.0, by = 0.5)
ylims2 <- c(0, 2.2)
y_title <- "Additional quality adjusted survival, \nQALYs"
p_qaly_prim <- .get_p_ly_undertreat(df = filter(gamma, secondary == 0), 
                                  title = title_no_cv, 
                                  xvar = xvar, yvar1 = yvar1, yvar2 = yvar2, 
                                  adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title, 
                                  bottom_graph = TRUE, left_graph = TRUE)
p_qaly_sec <- .get_p_ly_undertreat(df = filter(gamma, secondary == 1), 
                                 title = title_cv, 
                                 xvar = xvar, yvar1 = yvar1, yvar2 = yvar2, 
                                 adj = adj, ylims = ylims, ylims2 = ylims2, y_title = y_title, 
                                 bottom_graph = TRUE, right_graph = TRUE)

# combine
setwd(model_paper_dir)
png("ly_qaly_inc.png", height = 150+130, width = 220+200, units = "mm", res = 1000) 
grid.newpage()
pushViewport(viewport(
  height = unit(150+130, "mm"), width = unit(220+200, "mm"), 
  layout = grid.layout(2, 2, heights = unit(c(150, 130), "mm"), widths = unit(c(220, 200), "mm"))))
print(p_ly_prim, 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_ly_sec, 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_qaly_prim, 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_qaly_sec, 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()