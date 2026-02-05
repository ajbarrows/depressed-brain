# ---
# title: "multiverse analysis - repronim 2026"
# output:
#   html_document:
#     toc: true
#     toc_depth: 3
# ---

proj_root = "/Users/clarefmccann/University of Oregon Dropbox/Clare McCann/mine/projects/abcd-projs/PATHS/depressed-brain/"

pacman::p_load("dplyr")
source(paste0(proj_root, "R/PROCESS_v4/process.R"))

harm_data <- read.csv(paste0(proj_root, "data/processed/harmonized_tabular.csv")) %>% 
  mutate(age_centered = ab_g_dyn__visit_age - mean(ab_g_dyn__visit_age),
         sex_coded = ifelse(ab_g_stc__cohort_sex == "female", 1,
                            ifelse(ab_g_stc__cohort_sex == "male", 0, 
                                   ab_g_stc__cohort_sex)))
         #sex_coded = as.factor(sex_coded))
nonharm_data <- read.csv(paste0(proj_root, "data/processed/processed_tabular.csv")) %>% 
  mutate(age_centered = ab_g_dyn__visit_age - mean(ab_g_dyn__visit_age),
         sex_coded = ifelse(ab_g_stc__cohort_sex == "female", 1,
                            ifelse(ab_g_stc__cohort_sex == "male", 0, 
                                   ab_g_stc__cohort_sex)))

str(harm_data)
str(nonharm_data)

miss_by_wave <- harm_data %>%
  group_by(session_id) %>%
  summarise(across(
    where(~ TRUE),
    ~ mean(is.na(.)),
    .names = "miss_{.col}"
  ))

miss_by_wave

wave_counts <- harm_data %>%
  group_by(participant_id) %>%
  summarise(
    n_waves = n_distinct(session_id),
    .groups = "drop"
  )

n_all3 <- wave_counts %>%
  filter(n_waves == 3) %>%
  nrow()

n_all3

bl_harm_data <- harm_data %>%
  filter(session_id == "baseline") %>%
  select(participant_id,
         sex_coded,
         starts_with(c("baseline_maternal_", "age"))) %>% 
  rename(age_centered_bl = age_centered)

colMeans(is.na(bl_harm_data)) %>%  
  sort(decreasing = TRUE) %>%  
  head(20)

y2_harm_data <- harm_data %>%
  filter(session_id == "year_2") %>%
  select(participant_id,
         matches("^mr_y_.*\\.combat$"))

colMeans(is.na(y2_harm_data)) %>%  
  sort(decreasing = TRUE) %>%  
  head(20)

y4_harm_data <- harm_data  %>%
  filter(session_id == "year_4") %>%
  select(participant_id,
         matches("^mh_t_(bpm|total_ksads)"))

colMeans(is.na(y4_harm_data)) %>%  
  sort(decreasing = TRUE) %>%  
  head(20)

med_df_harm <- bl_harm_data %>%
  left_join(y2_harm_data, by = "participant_id") %>%
  left_join(y4_harm_data, by = "participant_id")

colMeans(is.na(med_df_harm)) %>%  
  sort(decreasing = TRUE) %>%  
  head(20)

vars <- c(
  "mh_t_bpm__int_sum",
  "baseline_maternal_asr_int",
  "mr_y_smri__vol__aseg__ag__lh_sum.combat",
  "age_centered_bl",
  "sex_coded"
)

sapply(med_df_harm[vars], class)

med_df_harm2 <- med_df_harm %>% 
  mutate(mh_t_bpm__int_sum = as.numeric(mh_t_bpm__int_sum),
         baseline_maternal_asr_int = as.numeric(baseline_maternal_asr_int),
         mr_y_smri__vol__aseg__ag__lh_sum.combat = as.numeric(mr_y_smri__vol__aseg__ag__lh_sum.combat),
         sex_coded = as.numeric(sex_coded),
         across(c(mh_t_bpm__int_sum, baseline_maternal_asr_int, mr_y_smri__vol__aseg__ag__lh_sum.combat), ~ as.numeric(scale(.x))))


################ HARMONIZED LEFT AMYGDALA W/O ICV ######################
#process(data = med_df_harm2, cov = "sex_coded", y = "mh_y_bpm__int_sum", x = "baseline_maternal_asr_int", m = "mr_y_smri__vol__aseg__ag__lh_sum.combat", boot = 10000, model = 4, total = 1)
library(lavaan)
model <- '
  mr_y_smri__vol__aseg__ag__lh_sum.combat ~ a*baseline_maternal_asr_int + s1*sex_coded
  mh_t_bpm__int_sum ~ b*mr_y_smri__vol__aseg__ag__lh_sum.combat + cprime*baseline_maternal_asr_int + s2*sex_coded

  indirect := a*b
  direct   := cprime
  total    := cprime + (a*b)
'

fit_fiml <- sem(
  model,
  data = med_df_harm2,
  missing = "fiml",     # the point of this exercise
  estimator = "MLR"     # robust SEs (nice default)
)

varTable(fit_fiml) %>%
  dplyr::arrange(desc(ov.var)) %>%
  dplyr::select(name, ov.var, nobs)


summary(fit_fiml, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
parameterEstimates(fit_fiml, standardized = TRUE) %>%
  dplyr::filter(label %in% c("a","b","cprime") | grepl("indirect|direct|total", lhs))