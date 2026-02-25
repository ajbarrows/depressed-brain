library(neuroCombat)

# CFM paths
#proj_path <- here::here()
#processed_path <- paste0(proj_path, "/data/processed/")

processed_path <- "./data/processed/"


load_subcortical_data <- function(fpath) {
  arrow::read_parquet(fpath) |>
    tidyr::drop_na(starts_with("mr_y_smri"))
}

pivot_to_matrix <- function(df) {
    df |>
        dplyr::select(participant_id, starts_with("mr_y_smri")) |>
        dplyr::distinct(participant_id, .keep_all = TRUE) |>
        tidyr::pivot_longer(-participant_id, names_to = 'feature') |>
        tidyr::pivot_wider(
            id_cols = feature,
            names_from = participant_id,
            values_from = value
        ) |>
        tibble::column_to_rownames(var = "feature")
}

pivot_to_dataframe <- function(combat_result) {
  combat_result$dat.combat |>
    as.data.frame() |>
    tibble::rownames_to_column() |>
    tidyr::pivot_longer(-rowname, names_to = 'participant_id') |>
    tidyr::pivot_wider(
      id_cols = participant_id,
      names_from = rowname,
      values_from = value
    )
}



harmonize_site <- function(df) {
  
  sub <- pivot_to_matrix(df)
  
  mod <- model.matrix(~ age + sex, data = df)
  batch <- as.numeric(as.factor(df$scanner_serial))

  subvolume_combat <- neuroCombat(dat = sub, batch = batch)
  
  # extract harmonized data
  harmonized_subcortical_volume <- pivot_to_dataframe(subvolume_combat)

  brainname_subvolume <- df |>
    dplyr::select(starts_with('mr_y_smri')) |>
    colnames()

  dplyr::left_join(
    harmonized_subcortical_volume,
    df |> dplyr::select(-all_of(brainname_subvolume)),
    by = 'participant_id'
  )
}




