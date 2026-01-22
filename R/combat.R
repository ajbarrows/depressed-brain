library(longCombat)

processed_path <- "./data/processed/"


load_subcortical_data <- function(fpath) {
  arrow::read_parquet(fpath) |>
    tidyr::drop_na()
}


harmonize_site <- function(df) {
  # brainname_subvolume <- c(names(df)[6:9])
  brainname_subvolume <- df |> dplyr::select(starts_with("mr_y_smri")) |> names()

  subvolume_combat <- longCombat(
    idvar = 'participant_id',
    timevar = 'session_id',
    # could also be smri_visitid
    batchvar = 'mr_y_adm__info__dev_serial',
    # should be a factor
    features = brainname_subvolume,
    formula = 'ab_g_dyn__visit_age + ab_g_stc__cohort_sex + session_id',
    ranef = '(1|participant_id)',
    data = df
  )


  # extract harmonized data
  harmonized_subcortical_volume <- subvolume_combat$data_combat

  dplyr::left_join(
    harmonized_subcortical_volume,
    df |> dplyr::select(-all_of(brainname_subvolume)),
    by = c(
      'participant_id',
      'session_id',
      'mr_y_adm__info__dev_serial'
    )
  )

}




