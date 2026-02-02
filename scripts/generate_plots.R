source("R/plot.R")

data_path <- "./data/processed/harmonized_tabular.parquet"
plot_path <- "./reports/figures/"

df <- arrow::read_parquet(data_path)

reshaped <- reshape_for_plotting(df)
brain_reshaped <- reshaped %>% reshape_for_brain_plotting()

plot_cbcl_by_sex(reshaped)
ggsave(paste0(plot_path, "mh_by_sex.png"), dpi = 'retina', height=4)

plot_mh_by_brain(brain_reshaped, structure_string = "Hippocampus")
ggsave(paste0(plot_path, "mh_by_hippocampus.png"), dpi = "retina", width=10, height=12)

plot_mh_by_brain(brain_reshaped, structure_string = "Amygdala")
ggsave(paste0(plot_path, "mh_by_amygdala.png"), dpi = "retina", width=10, height=12)

make_pairplot(df, outcomes_map, "Depression Symptoms (all time points)")
ggsave("./reports/figures/mh_correlation_all.png", dpi = 'retina', width=8, height=8)

tpts <- c("baseline", "year_2", "year_4", "year_6")
fpath <- "./reports/figures/"
for (tpt in tpts) {
  title <- paste0("Depression Symptoms (", tpt, ")")
  fname <- paste0(fpath, "mh_correlation_", tpt, ".png")
    
  p <- df %>%
    filter(session_id == tpt) %>%
    make_pairplot(outcomes_map, title)

  ggsave(fname, plot = p, dpi = "retina", width=8, height=8)
}