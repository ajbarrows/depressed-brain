source("R/plot.R")

data_path <- "./data/processed/harmonized_tabular.parquet"
plot_path <- "./reports/figures/"

df <- arrow::read_parquet(data_path)

reshaped <- reshape_for_plotting(df)
brain_reshaped <- reshaped %>% reshape_for_brain_plotting()

plot_cbcl_by_sex(reshaped)
ggsave(paste0(plot_path, "cbcl_by_sex.png"), dpi = 'retina', height=4)

plot_mh_by_brain(brain_reshaped, structure_string = "Hippocampus")
ggsave(paste0(plot_path, "mh_by_hippocampus.png"), dpi = "retina", width=10)

plot_mh_by_brain(brain_reshaped, structure_string = "Amygdala")
ggsave(paste0(plot_path, "mh_by_amygdala.png"), dpi = "retina", width=10)
