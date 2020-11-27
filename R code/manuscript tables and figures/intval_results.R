## POINT ESTIMATES -------------------------------------------

rm(list = ls())

library(here)

library(data.table)
library(tidyverse)

### directories and functions

source(here("src", "config", "define_dir.R"))
source(here("src", "model_paper", "model_paper_functions.R"))

### load baseline data ########

setwd(data_processed_dir)
df_baseline_all <- get(load("df_baseline_imputed.Rdata"))
df_baseline <- select(df_baseline_all, patid, ckd) %>%
  mutate(ckd2 = ifelse(ckd %in% c("ckd4", "ckd5"), "ckd4/5", ckd)) %>%
  mutate(ckd2 = factor(ckd2, labels=c("ckd2", "ckd3a", "ckd3b", "ckd4/5"))) %>%
  select(-ckd) %>%
  rename(ckd = ckd2)

### prepare observed data ########

setwd(data_processed_dir)
df_obs_all <- get(load("df_Aevents_Cross_Validation.Rdata"))

# leave estimation cohort only
df_obs <- filter(df_obs_all, estimation == 0)
df_id_all <- get(load("df_id_all.Rdata"))
df_obs <- merge(df_obs, select(df_id_all, patid, patid_temp))
id <- unique(df_obs$patid)
id_temp <- unique(df_obs$patid_temp)

# clean
df_obs <- df_obs %>%
  mutate(TTE_HF_hosp2 = ifelse(HF_hosp == 0 & prev_HF == 1, 0, TTE_HF_hosp)) %>%
  rename(endpt_nvd = NVD, fu_nvd = TTE_NVD, endpt_vd = VD, fu_vd = TTE_VD,
         endpt_vd_stroke = Stroke_or_VD, fu_vd_stroke = TTE_Stroke_or_VD,
         endpt_vd_stroke_mi = MI_Stroke_or_VD, fu_vd_stroke_mi = TTE_MI_Stroke_or_VD,
         endpt_hf = HF_hosp, fu_hf = TTE_HF_hosp2) %>%
  select(patid, patid_temp, cycle, starts_with("endpt_"), starts_with("fu_"))

df_obs <- merge(df_obs, df_baseline, by = "patid")

cols_to_summarise <- setdiff(colnames(df_obs), c("patid", "patid_temp", "cycle", "ckd"))

df_obs <- df_obs %>%
  select(-patid) %>%
  arrange(patid_temp) %>%
  group_by(patid_temp, ckd) %>%
  summarise_at(cols_to_summarise, sum) %>%
  ungroup()

### prepare predicted data ########

df_pred_all <- .combine_4_into_1(file_prefix = "df_output_int_val", dir = datasets_output_dir)
df_pred <- df_pred_all %>%
  filter(patid_temp %in% id_temp) %>%
  rename(endpt_nvd = endpt_NVD, endpt_vd = endpt_VD, endpt_hf = HF_hosp) %>%
  mutate(fu_nvd = alive_at_start, fu_vd = alive_at_start - endpt_nvd) %>% 
  mutate(fu_vd_stroke = shift(fu_s_vd_temp, 1), fu_vd_stroke_mi = shift(fu_mi_s_vd_temp, 1)) %>% 
  mutate(fu_vd_stroke = ifelse(cycle == 1, 1, fu_vd_stroke), 
         fu_vd_stroke_mi = ifelse(is.na(fu_vd_stroke_mi), 1, fu_vd_stroke_mi)) %>%
  mutate(fu_vd_stroke = fu_vd_stroke - endpt_nvd, 
         fu_vd_stroke_mi = fu_vd_stroke_mi - endpt_nvd) %>%
  mutate(fu_hf = alive_at_start - prev_HF) %>%
  select(-ends_with("vd_temp")) %>%
  select(patid_temp, cycle, starts_with("endpt_"), starts_with("fu_")) %>% 
  arrange(patid_temp)

df_pred <- merge(df_pred, select(df_id_all, patid, patid_temp), by = "patid_temp")
df_pred <- merge(df_pred, df_baseline, by = "patid")
df_pred <- df_pred %>%
  select(-patid) %>%
  arrange(patid_temp, cycle)

### internal validation ##################

library(survival)
library(ggplot2)
library(scales)
library(graphics)
library(grDevices)

### calculate cumulative rates

endpts <- list(
  list(lab = "Non-vascular \ndeath", endptAct = "endpt_nvd", endptPred = "endpt_nvd", fuvarAct="fu_nvd", fuvarPred="fu_nvd"),
  list(lab = "Vascular \ndeath", endptAct = "endpt_vd", endptPred = "endpt_vd", fuvarAct="fu_vd", fuvarPred="fu_vd"),
  list(lab = "VD or Stroke", endptAct = "endpt_vd_stroke", endptPred = "endpt_vd_stroke", fuvarAct="fu_vd_stroke", fuvarPred="fu_vd_stroke"),
  list(lab = "MI,VD \nor Stroke", endptAct = "endpt_vd_stroke_mi", endptPred = "endpt_vd_stroke_mi", fuvarAct="fu_vd_stroke_mi", fuvarPred="fu_vd_stroke_mi"),
  list(lab = "Heart Failure \nHospitalisation", endptAct = "endpt_hf", endptPred = "endpt_hf", fuvarAct="fu_hf", fuvarPred="fu_hf"))

subgroups <- levels(df_pred$ckd)
Y <- 5
output_act <- .cumprob_Act(dfAct = df_obs, Y = Y, subgroups = subgroups, endpts = endpts)
output_pred <- .cumprob_Pred(dfPred = df_pred, subgroups = subgroups, endpts = endpts, Y = Y)
output <- merge(output_act, output_pred, by=c("subgroup", "lab", "year"))

endpt <- endpts[[1]]
subgroup <- subgroups[1]

# graphical parameters

subgroup2_levels <- c("G2", "G3a", "G3b", "G4/G5")
subgroup2_labels <- c(
  paste0("eGFR 60", "\U2013", "89 mL/min/1.73", "\U00B2", "\n\n(G2)"),
  paste0("eGFR 45", "\U2013", "59 mL/min/1.73", "\U00B2", "\n\n(G3a)"),
  paste0("eGFR 30", "\U2013", "44 mL/min/1.73", "\U00B2", "\n\n(G3b)"),
  paste0("eGFR \u2264 29mL/min/1.73", "\U00B2", "\nnot on RRT\n(G4/5)")
)

lab_levels <- c("Non-vascular \ndeath", "Vascular \ndeath", "VD or Stroke", "MI,VD \nor Stroke", "Heart Failure \nHospitalisation")
lab_labels <- c("Nonvascular \ndeath", 
                "Vascular death", 
                "Vascular death \nor\nstroke", 
                "Vascular death,\nstroke\nor\nmyocardial infarction", 
                "Heart failure \nhospitalisation")
output <- mutate(output, 
                 lab = factor(lab, levels = lab_levels, labels = lab_labels),
                 subgroup2 = factor(subgroup, labels = subgroup2_labels))

p_breaks <- seq(from = 0, to = 0.3, by = 0.1)
x_axis_label <- "year of follow-up in CPRD validation cohort"
y_axis_label <- ""

setwd(model_paper_dir)

### cardiovascular endpoints (main paper)

df <- filter(output, lab !="Nonvascular \ndeath" & lab != "Heart failure \nhospitalisation")
p <- .get_p_intval(df = df, year = year, p_breaks = p_breaks, x_axis_label = x_axis_label, y_axis_label = y_axis_label)

p

ggsave(p, filename="internal_validation_MAIN.png", height = 245, width = 320, units = "mm", dpi = 300)

# NVD and HF (supplementary)

df <- filter(output, lab %in% c("Nonvascular \ndeath", "Heart failure \nhospitalisation"))
p <- .get_p_intval(df = df, year = year, p_breaks = p_breaks, x_axis_label = x_axis_label, y_axis_label = y_axis_label)
ggsave(p, filename="internal_validation_SUPPL.png", height = 150, width = 320, units = "mm", dpi = 300)

### validation by region ################

subgroup2_labels <- c(
  paste0("eGFR 60", "\U2013", "89", "\n mL/min/1.73", "\U00B2", "\n\n(G2)"),
  paste0("eGFR 45", "\U2013", "59", "\n mL/min/1.73", "\U00B2", "\n\n(G3a)"),
  paste0("eGFR 30", "\U2013", "44", "\n mL/min/1.73", "\U00B2", "\n\n(G3b)"),
  paste0("eGFR \u2264 29", "\nmL/min/1.73", "\U00B2", "\nnot on RRT\n(G4/5)")
)


# define correct datasets

# add region data 
df_practice <- get(load("K:\\SHARP_EE\\Monitoring CKD CHF\\Data\\Derived CPRD data\\df_practice.Rdata"))
df_loc <- merge(select(df_baseline_all, patid, pracid), select(df_practice, pracid, region_txt), by = "pracid", all.x = TRUE)
nrow(df_loc)
df_loc <- merge(df_loc, select(df_id_all, patid, patid_temp))
df_loc <- df_loc %>%
  mutate(region = region_txt) %>%
  select(-patid) %>%
  select(patid_temp, region)
table(df_loc$region, exclude = NULL)
df_obs <- merge(df_obs, df_loc)
df_obs <- arrange(df_obs, patid_temp)
df_pred <- merge(df_pred, df_loc)
df_pred <- arrange(df_pred, patid_temp, cycle)

# endpoints

lab_levels <- c("MI,VD \nor Stroke")
lab_labels <- c("Vascular death,\nstroke\nor\nmyocardial infarction")

# produce a graph for each region

Y <- 5
p_breaks <- seq(from = 0, to = 0.3, by = 0.1)
x_axis_label <- "year of follow-up in CPRD validation cohort"
y_axis_label <- ""

regions <- sort(unique(df_loc$region))

p_list <- list()

for (region_temp in regions) {
  
  # filter the datasets
  df_obs_temp <- filter(df_obs, region == region_temp)
  df_pred_temp <- filter(df_pred, region == region_temp)
  
  # summarise the results
  output_obs <- .cumprob_Act(dfAct = df_obs_temp, subgroups = subgroups, endpts = endpts, Y = Y)
  # predicted
  output_pred <- .cumprob_Pred(dfPred = df_pred_temp, subgroups = subgroups, endpts = endpts, Y = Y)
  # merge
  output <- merge(output_obs, output_pred, by=c("subgroup", "lab", "year"))
  output <- mutate(output, 
                   lab = factor(lab, levels = lab_levels, labels = lab_labels),
                   subgroup2 = factor(subgroup, labels = subgroup2_labels))
  df <- filter(output, lab == "Vascular death,\nstroke\nor\nmyocardial infarction")
  p <- .get_p_intval(df = df, year = year, p_breaks = p_breaks, x_axis_label = x_axis_label, y_axis_label = y_axis_label)
  p <- p + ggtitle(region_temp) + theme(strip.text.y = element_blank(),
                                        #axis.title.x = element_blank(),
                                        plot.title = element_text(hjust = 0.5, size = 25))
  p_list[[region_temp]] <- p
}

library(grid)
library(gridExtra)

# combine
setwd(model_paper_dir)
png("validation_region_north.png", height = 200,  width = 600, units = "mm", res = 300)
grid.newpage()
pushViewport(viewport(
  height = unit(200, "mm"), width = unit(600, "mm"), 
  layout = grid.layout(1, 3, heights = unit(rep(200, 1), "mm"), width = unit(rep(200, 3), "mm"))))

print(p_list[["North West"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_list[["North East"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_list[["Yorkshire & The Humber"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 3))

dev.off()

png("validation_region_mid.png", height = 200,  width = 600, units = "mm", res = 300)
grid.newpage()
grid.newpage()
pushViewport(viewport(
  height = unit(200, "mm"), width = unit(600, "mm"), 
  layout = grid.layout(1, 3, heights = unit(rep(200, 1), "mm"), width = unit(rep(200, 3), "mm"))))

print(p_list[["West Midlands"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_list[["East Midlands"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_list[["East of England"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 3))

dev.off()

png("validation_region_south.png", height = 400,  width = 400, units = "mm", res = 300)
grid.newpage()
pushViewport(viewport(
  height = unit(400, "mm"), width = unit(400, "mm"), 
  layout = grid.layout(2, 2, heights = unit(rep(200, 2), "mm"), width = unit(rep(200, 2), "mm"))))

print(p_list[["London"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_list[["South West"]], 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_list[["South Central"]], 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p_list[["South East Coast"]], 
      vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

dev.off()
