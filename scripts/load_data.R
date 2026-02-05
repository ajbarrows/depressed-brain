source("/Users/clarefmccann/University of Oregon Dropbox/Clare McCann/mine/projects/abcd-projs/PATHS/depressed-brain/R/data.R")

data <- load_data(
  repeated_tables,
  repeated_key,
  single_tpt_tables,
  single_tpt_key,
  timepoint_map,
  sex_map
)

data$df %>%
  arrow::write_parquet(paste0(processed_path, "full_dataset.parquet"))

data$subset %>%
  arrow::write_parquet(paste0(processed_path, "processed_tabular.parquet"))
data$subset %>%
  write.csv(paste0(processed_path, "processed_tabular.csv"))
