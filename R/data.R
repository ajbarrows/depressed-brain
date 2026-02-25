library(dplyr)

## CFM paths
#proj_path <- here::here()
#data_root <- Sys.getenv("HOME_DIR")

#data_path <- paste0(data_root, "/projects/abcd-projs/abcd-data-release-6.0/nbdc-tools-data/")
#processed_path <- paste0(proj_path, "/data/processed/")

data_path <- "./data/raw/"
processed_path <- "./data/processed/"

repeated_tables <- list(
  "mr_y_smri__vol__aseg" = c(
    "mr_y_smri__vol__aseg__ag__lh_sum",
    "mr_y_smri__vol__aseg__ag__rh_sum",
    "mr_y_smri__vol__aseg__hc__lh_sum",
    "mr_y_smri__vol__aseg__hc__rh_sum",
    "mr_y_smri__vol__aseg__icv_sum"
  ),
  "mr_y_adm__info" = c("mr_y_adm__info__dev_serial"),
  "mr_y_qc__incl" = c("mr_y_qc__incl__smri__t1_indicator"),
  "le_l_prenatal" = c(
    "le_l_prenatal__addr1__no2_mean",
    "le_l_prenatal__addr1__pm25_mean",
    "le_l_prenatal__addr1__o3_mean"
  ),

  "mh_p_abcl" = c("mh_p_abcl__cg2_001"),
  "mh_p_asr" = c(
    "mh_p_asr_dtt",
    "mh_p_asr__synd__int_sum",
    "mh_p_asr__synd__anxdep_sum"
  ),
  'ab_g_dyn' = NULL,
  "mh_p_cbcl" = c(
    "mh_p_cbcl__synd__attn_sum",
    "mh_p_cbcl__synd__anxdep_sum",
    "mh_p_cbcl__synd__int_sum",
    "mh_p_cbcl__synd__wthdep_sum",
    "mh_p_cbcl__dsm__dep_sum"
  ),
  'mh_y_ksads__dep' = NULL,
  'mh_t_bpm' = c('mh_t_bpm__int_sum'),
  'mh_y_bpm' = c('mh_y_bpm__int_sum')
)

single_tpt_tables <- list(
  "ab_g_stc" = c(
    "ab_g_stc__cohort_sex",
    "ab_g_stc__design_id__fam"
  )
)

repeated_key <- c("participant_id", "session_id")
single_tpt_key <- "participant_id"

timepoint_map <- c(
  "ses-00A" = "baseline",
  "ses-02A" = "year_2",
  "ses-04A" = "year_4",
  "ses-06A" = "year_6"
)

sex_map <- c(
  "2" = "female",
  "1" = "male"
)

load_and_join <- function(variables, key) {
  tables <- purrr::imap(
    variables,
    ~ {
      fpath <- paste0(data_path, .y, ".parquet")
      if (is.null(.x)) {
        arrow::read_parquet(fpath)
      } else {
        arrow::read_parquet(fpath, col_select = c(all_of(key), .x))
      }
    }
  )

  purrr::reduce(tables, dplyr::left_join, by = key)
}

filter_mri_qc <- function(df) {
  df %>%
    filter(mr_y_qc__incl__smri__t1_indicator == 1)
}

# figure out what day the ASR was collected
extract_corresponding_day <- function(df) {
  df %>%
    mutate(
      ab_g_dyn__visit__day1_dt = ab_g_dyn__visit_dtt
    ) %>%
    select(
      participant_id,
      session_id,
      contains("dt"),
      -ab_g_dyn__visit_dtt
    ) %>%
    rowwise() %>%
    mutate(
      num_visit_days = sum(!is.na(c_across(ab_g_dyn__visit__day1_dt:ab_g_dyn__visit__day4_dt)))
    ) %>%
    ungroup() %>%
    mutate(across(where(lubridate::is.POSIXct), lubridate::date)) %>%
    tidyr::pivot_longer(
      cols = ends_with("dt"),
      values_to = 'date',
      names_to = 'date_var'
    ) %>%
    filter(mh_p_asr_dtt == date | num_visit_days  == 1) %>%
    mutate(asr_day = stringr::str_extract(date_var, "day\\d+"))
}
limit_to_bio_mom <- function(df, timepoint = 'ses-00A') {

  asr_admin_date <- extract_corresponding_day(df)

  bio_mom <- df %>%
    select(participant_id, session_id, ends_with("inform")) %>%
    tidyr::pivot_longer(
      cols = ends_with('inform')
    )  %>%
    mutate(inform_day = stringr::str_extract(name, "day\\d+")) %>%
    left_join(
      asr_admin_date,
      by = c("participant_id", "session_id"),
      relationship = "many-to-many"
    )  %>%
    filter(inform_day == asr_day | num_visit_days  == 1) %>%
    distinct(participant_id, session_id, .keep_all = TRUE) %>%
    select(
      participant_id,
      session_id,
      informant = value
    ) %>%
    filter(informant == 1) %>% # biological mother
    filter(session_id == timepoint)
  # #

  #
  sub <- df %>%
    right_join(bio_mom, by = c("participant_id", "session_id")) %>%
    filter(session_id == timepoint) %>%
    select(
      participant_id,
      baseline_maternal_asr_int = mh_p_asr__synd__int_sum,
      baseline_maternal_asr_anxdep = mh_p_asr__synd__anxdep_sum
    )

  df %>%
    right_join(sub, by = "participant_id")

}

summarize_ksads_dep_sx <- function(df) {
  df %>%
    select(
      (starts_with('mh_y_ksads') & ends_with('sx') & contains("pres")) # grab symptom scales (present only)
    ) %>%
      mutate(across(everything(), ~case_match(
        .x,
        "1" ~ "1",
        "0" ~ "0",
        .default = NA_character_
      ))) %>%
      mutate(across(everything(), as.numeric)) %>%
      mutate(total_ksads_dep_sx = rowSums(across(everything()), na.rm=TRUE)) %>%
      select(total_ksads_dep_sx) %>%
      bind_cols(
        df %>%
          select(
            participant_id, 
            session_id,
            mh_y_ksads__dep__mdd__pres_dx,
            mh_y_ksads__dep__pdd__pres_dx,
            !starts_with('mh_y_ksads')
          )
      ) %>%
      select(participant_id, session_id, everything())
}


# format_subset <- function(df, timepoint_map, sex_map) {
#   df %>%
#     select(
#       participant_id,
#       session_id,
#       # sex = ab_g_stc__cohort_sex,
#       # age = ab_g_dyn__visit_age,
#       # site = ab_g_dyn__design_site,
#       # scanner_serial = ab_g_dyn__design_mr__serial,
#       # family_id = ab_g_stc__design_id__fam,
#       # anxdep = mh_p_cbcl__synd__anxdep_sum,
#       # internal = mh_p_cbcl__synd__int_sum,
#       # attention = mh_p_cbcl__synd__attn_sum,
#       # left_amygdala_vol = mr_y_smri__vol__aseg__ag__lh_sum,
#       # right_amygdala_vol = mr_y_smri__vol__aseg__ag__rh_sum,
#       # left_hippocampus_vol = mr_y_smri__vol__aseg__hc__lh_sum,
#       # right_hippocampus_vol = mr_y_smri__vol__aseg__hc__rh_sum,
#       ab_g_stc__cohort_sex,
#       ab_g_dyn__visit_age,
#       ab_g_dyn__design_site,
#       ab_g_stc__design_id__fam,
#       mh_p_cbcl__synd__anxdep_sum,
#       mh_p_cbcl__synd__int_sum,
#       mh_p_cbcl__synd__attn_sum,
#       mh_p_cbcl__synd__wthdep_sum,
#       mh_p_cbcl__dsm__dep_sum,
#       mh_t_bpm__int_sum,
#       mh_y_bpm__int_sum,
#       mr_y_smri__vol__aseg__ag__lh_sum,
#       mr_y_smri__vol__aseg__ag__rh_sum,
#       mr_y_smri__vol__aseg__hc__lh_sum,
#       mr_y_smri__vol__aseg__hc__rh_sum,
#       mr_y_smri__vol__aseg__icv_sum,
#       mr_y_adm__info__dev_serial,
#       baseline_maternal_asr_int,
#       baseline_maternal_asr_anxdep,
#       mh_y_ksads__dep__mdd__pres_dx,
#       mh_y_ksads__dep__pdd__pres_dx,
#       total_ksads_dep_sx
#     ) %>%
#     mutate(
#       session_id = recode(session_id, !!!timepoint_map),
#       ab_g_stc__cohort_sex = recode(ab_g_stc__cohort_sex, !!!sex_map),
#       mr_y_adm__info__dev_serial = factor(mr_y_adm__info__dev_serial)
#     ) %>%
#     arrange(participant_id, session_id) %>%
#     tidyr::drop_na(!(starts_with("mh_y_bpm")))
# }


subset_baseline <- function(df) {
  df %>%
    filter(session_id == 'baseline') %>%
    select(
      participant_id,
      sex,
      age = ab_g_dyn__visit_age,
      site = ab_g_dyn__design_site,
      scanner_serial,
      family_id = ab_g_stc__design_id__fam,
      baseline_maternal_asr_int,
      baseline_maternal_asr_anxdep
    )
}

subset_timepoint <- function(
  df,
  timepoints = c("year_2"),
  var_string = "mr_y_smri"
) {
  n_vars <- df %>% select(contains(var_string)) %>% ncol()

  tmp <- df %>%
    filter(session_id %in% timepoints) %>%
    filter_mri_qc() %>%
    select(
      participant_id,
      session_id,
      contains(var_string)
    ) %>%
    tidyr::pivot_wider(
      id_cols = participant_id,
      names_from = session_id,
      values_from = contains(var_string)
    )

  if (n_vars < 2) {
    return(
      tmp %>%
        rename_with(
          ~ paste0(var_string, "_", .x),
          .cols = -c(participant_id)
        )
    )
  } else {
    return(tmp)
  }
}

format_subset <- function(df, timepoint_map, sex_map) {
  tmp <- df %>%
    mutate(
      session_id = recode(session_id, !!!timepoint_map),
      sex = recode(ab_g_stc__cohort_sex, !!!sex_map),
      scanner_serial = factor(mr_y_adm__info__dev_serial)
    )
  tmp %>%
    subset_baseline() %>%
    join_tables(
      subset_timepoint(
        tmp,
        timepoints = c("year_2"),
        var_string = "mr_y_smri"
      )
    ) %>%
    join_tables(
      subset_timepoint(
        tmp,
        timepoints = c("year_4"),
        var_string = "mh_y_bpm__int"
      )
    ) %>%
      tidyr::drop_na()
}


join_tables <- function(df, other, key = 'participant_id') {
  df %>%
    left_join(other, by = key)
}


load_data <- function(
    repeated_tables,
    repeated_key,
    single_tpt_tables,
    single_tpt_key,
    timepoint_map,
    sex_map
) {

  df <- load_and_join(repeated_tables, repeated_key) %>%
    join_tables(load_and_join(single_tpt_tables, single_tpt_key))

  subset <- df %>%
    filter_mri_qc() %>%
    limit_to_bio_mom() %>%
    summarize_ksads_dep_sx() %>%
    format_subset(timepoint_map, sex_map)

  return(list(
    "df" = df,
    "subset" = subset
  ))
}
