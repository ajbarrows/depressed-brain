library(ggplot2)
library(dplyr)

outcomes_map <- c(
  "mh_p_cbcl__synd__wthdep_sum" = "CBCL Withdrawn/Depressed",
  "mh_p_cbcl__synd__int_sum" = "CBCL Internalizing",
  'mh_y_bpm__int_sum' = "BPM Internalizing (Youth)",
  'total_ksads_dep_sx' = "KSADS Dep. Symptoms (sum)"
)

brain_map <- c(
  "mr_y_smri__vol__aseg__ag__lh_sum.combat" = "Left Amygdala Volume",
  "mr_y_smri__vol__aseg__ag__rh_sum.combat" = "Right Amygdala Volume",
  "mr_y_smri__vol__aseg__hc__lh_sum.combat" = "Left Hippocampus Volume",
  "mr_y_smri__vol__aseg__hc__rh_sum.combat" = "Right Hippocampus Volume"
)

##### manipulation functions #####

join_median_timepoint_age <- function(df) {
  summary <- df %>%
    group_by(session_id) %>%
    summarize(median_age = median(ab_g_dyn__visit_age))
  
  df %>%
    left_join(summary, by = "session_id")
}

mask_outliers <- function(df, group_col, value_col, masked_col = paste0(value_col, "_masked")) {
  stats <- df %>%
    group_by(.data[[group_col]]) %>%
    summarize(
      .mean = mean(.data[[value_col]], na.rm = TRUE),
      .twosd = sd(.data[[value_col]], na.rm = TRUE) * 2,
      .groups = "drop"
    )

  df %>%
    left_join(stats, by = group_col) %>%
    mutate(!!masked_col := ifelse(
      .data[[value_col]] > .mean + .twosd |
        .data[[value_col]] <= .mean - .twosd,
      NA, .data[[value_col]])) %>%
    select(-c(.mean, .twosd))
}

reshape_for_plotting <- function(df) {
  df %>%
    join_median_timepoint_age() %>%
    tidyr::pivot_longer(
      cols = c(
        mh_p_cbcl__synd__wthdep_sum,
        mh_p_cbcl__synd__int_sum,
        mh_y_bpm__int_sum,
        total_ksads_dep_sx
      ),
      names_to = 'summary_score',
      values_to = 'summary_score_value'
    ) %>%
    mutate(summary_score = recode(summary_score, !!!(outcomes_map))) %>%
    mask_outliers("summary_score", "summary_score_value", "masked_summary_score") %>%
    select(
      age = ab_g_dyn__visit_age,
      sex = ab_g_stc__cohort_sex,
      everything()
    )
}

reshape_for_brain_plotting <- function(reshaped) {
  reshaped %>%
    select(
      participant_id,
      session_id,
      age,
      summary_score,
      summary_score_value,
      ends_with("combat"),
      -ends_with('serial'),
      -contains("icv")) %>%
    tidyr::pivot_longer(
      cols = ends_with("combat"),
      names_to = 'structure',
      values_to = 'volume'
    ) %>%
    mutate(structure_var = recode(structure, !!!(brain_map))) %>%
    mask_outliers("structure_var", "volume", "volume_masked")
}


##### plotting functions #####

plot_cbcl_by_sex <- function(plot_df) {
  ggplot(plot_df, aes(x = age, y = masked_summary_score)) +
    geom_point(alpha = 0.01) +
    geom_vline(
      aes(xintercept = median_age, color = factor(median_age)),
      linetype = "dashed") +
    facet_wrap(~summary_score, scales='free') +
    geom_smooth(aes(color = sex)) +
  scale_color_discrete(
    labels = setNames(
      c("Baseline", "Year 2", "Year 4", "Year 6"),
      sort(unique(plot_df$median_age))
    )
  ) +
    labs(
      color = "",
      x = "Age (years)",
      y = "Summary Scale Value",
      caption = "values >= 2sd above the mean removed for plotting"
    ) +
    theme_minimal()
}

plot_mh_by_brain <- function(
    df, 
    structure_string = 'Hippocampus', 
    alpha = 0.25, 
    transform = "log1p",
    color_label = "Log Summary Score Value"){
  df %>%
    filter(stringr::str_detect(structure_var, structure_string)) %>%
    ggplot(aes(x = age, y = volume_masked, color = summary_score_value)) +
      geom_point(alpha = alpha) +
      scale_color_viridis_c(transform = transform) +
      facet_grid(
        rows = vars(summary_score),
        cols = vars(structure_var)
      ) +
      theme_minimal(base_size = 15) +
      theme(
        legend.position = "top",
        legend.key.width = unit(2, "cm"),   # make legend bar wider
        legend.text = element_text(size = 10)  # bigger text
      ) +
      labs(
        color = color_label,
        y = "Volume",
        x = "Age (years)"
      ) 
}

make_pairplot <- function(df, outcomces_map, title = "") {
  scale_agreement <- function(data, mapping) {
    x_var <- rlang::as_name(mapping$x)
    y_var <- rlang::as_name(mapping$y)
    
    cor_val <- cor(data[[x_var]], data[[y_var]], use = "complete.obs")
    
    ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.05) +
      geom_smooth(method = 'lm') +
      annotate("text", x = Inf, y = Inf, 
              label = sprintf("r = %.2f", cor_val),
              hjust = 1.1, vjust = 1.5)
    }

  df %>%
    select(
      mh_y_bpm__int_sum,
      total_ksads_dep_sx,
      mh_p_cbcl__synd__wthdep_sum,
      mh_p_cbcl__synd__int_sum
    ) %>%
      rename_with(everything(), .fn = ~outcomes_map) %>%
      GGally::ggpairs(
        lower = list(continuous = scale_agreement),
        upper = list(continuous = "blank"),
        diag = list(continuous = "barDiag"),
        switch = "both"
      ) +
        theme_minimal(base_size = 13) +
        labs(title = title)
}