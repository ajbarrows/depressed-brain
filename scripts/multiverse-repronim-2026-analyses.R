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

harm_data <- read.csv(paste0(proj_root, "data/processed/harmonized_tabular.csv"))
nonharm_data <- read.csv(paste0(proj_root, "data/processed/processed_tabular.csv"))

process(data = harm_data, cov = c(), y = "", x = "", m = c(""), boot = 10000, model = 6, total = 1)
