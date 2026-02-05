source("/Users/clarefmccann/University of Oregon Dropbox/Clare McCann/mine/projects/abcd-projs/PATHS/depressed-brain/R/combat.R")

set.seed(42)

included_data <- load_subcortical_data(paste0(processed_path, "processed_tabular.parquet"))

harmonized <- harmonize_site(included_data)

harmonized |>
  arrow::write_parquet(paste0(processed_path, "harmonized_tabular.parquet"))

harmonized |>
  write.csv(paste0(processed_path, "harmonized_tabular.csv"))

