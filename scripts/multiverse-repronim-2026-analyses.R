# ---
# title: "multiverse analysis - repronim 2026"
# output:
#   html_document:
#     toc: true
#     toc_depth: 3
# ---

# --- project + deps ----------------------------------------------------

proj_root <- here::here()

source(file.path(proj_root, "R/PROCESS_v4/process.R"))

pacman::p_load(
  "arrow",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "htmltools",
  "stringr",
  "ggplot2",
  "forcats",
  "patchwork",
  "readr"
)

harm_data <- read_parquet(file.path(
  proj_root,
  "data/processed/harmonized_tabular.parquet"
)) %>%
  mutate(
    age_centered = age - mean(age),
    sex_coded = ifelse(sex == "female", 1, ifelse(sex == "male", 0, sex))
  )

nonharm_data <- read_parquet(file.path(
  proj_root,
  "data/processed/processed_tabular.parquet"
)) %>%
  mutate(
    age_centered = age - mean(age),
    sex_coded = ifelse(sex == "female", 1, ifelse(sex == "male", 0, sex))
  )

# --- model spec --------------------------------------------------------

data_map <- list(
  harmonized = harm_data,
  nonharmonized = nonharm_data
)

y_var <- "mh_y_bpm__int_year_4"
x_var <- "baseline_maternal_asr_int"
cov_base <- "sex_coded"
icv_var <- "mr_y_smri__vol__aseg__icv_sum_year_2"

spec <- tribble(
  ~structure    , ~side , ~m_var                                    ,
  "amygdala"    , "lh"  , "mr_y_smri__vol__aseg__ag__lh_sum_year_2" ,
  "amygdala"    , "rh"  , "mr_y_smri__vol__aseg__ag__rh_sum_year_2" ,
  "hippocampus" , "lh"  , "mr_y_smri__vol__aseg__hc__lh_sum_year_2" ,
  "hippocampus" , "rh"  , "mr_y_smri__vol__aseg__hc__rh_sum_year_2"
) |>
  crossing(
    dataset = names(data_map),
    icv = c(FALSE, TRUE)
  ) |>
  mutate(
    cov = map(icv, ~ if (.x) c(cov_base, icv_var) else cov_base),
    model_id = paste(
      dataset,
      structure,
      side,
      ifelse(icv, "with_icv", "no_icv"),
      sep = "_"
    )
  )

# --- output ------------------------------------------------------------

split_process_runs <- function(lines) {
  idx <- grep("^\\*{10,}\\s*PROCESS for R", lines)
  if (length(idx) == 0) {
    return(list(lines))
  }
  idx <- c(idx, length(lines) + 1)

  runs <- lapply(seq_len(length(idx) - 1), function(i) {
    lines[idx[i]:(idx[i + 1] - 1)]
  })
  runs[lengths(runs) > 0]
}

remove_boot_progress <- function(lines) {
  start <- grep("^Bootstrapping progress:", lines)
  if (length(start) == 0) {
    return(lines)
  }
  start <- start[1]

  after <- lines[start:length(lines)]
  end_rel <- grep("100%\\s*$", after)
  if (length(end_rel) == 0) {
    return(lines[-start])
  }

  end <- start + end_rel[1] - 1
  lines[-(start:end)]
}

extract_results_chunk <- function(lines) {
  start <- grep("^\\*{5,}\\s*Outcome Variable:", lines)
  if (length(start) == 0) {
    start <- grep("^Outcome Variable:", lines)
  }
  if (length(start) == 0) {
    return(character())
  }
  start <- start[1]

  end <- grep("^\\*{5,}\\s*ANALYSIS NOTES AND ERRORS", lines)
  end <- end[end > start]
  if (length(end) > 0) {
    end <- end[1] - 1
    chunk <- lines[start:end]
    chunk <- remove_boot_progress(chunk)
    return(chunk)
  }

  chunk <- lines[start:length(lines)]
  chunk <- remove_boot_progress(chunk)
  chunk
}

capture_process <- function(data, cov, y, x, m, boot = 10000) {
  all_lines <- capture.output(
    process(
      data = data,
      cov = cov,
      y = y,
      x = x,
      m = m,
      boot = boot,
      model = 4,
      total = 1
    )
  )

  runs <- split_process_runs(all_lines)

  cleaned <- lapply(runs, function(r) {
    chunk <- extract_results_chunk(r)
    chunk <- chunk[nzchar(chunk)]
    paste(chunk, collapse = "\n")
  })

  paste(cleaned, collapse = "\n\n-----\n\n")
}

# --- run models --------------------------------------------------------

models <- spec |>
  mutate(
    output = pmap_chr(
      list(
        data = map(dataset, ~ data_map[[.x]]),
        cov = cov,
        y = rep(y_var, n()),
        x = rep(x_var, n()),
        m = m_var
      ),
      capture_process
    )
  )

# --- parse process output ----------------------------------------------

parse_process_output <- function(output_text) {
  num_pat <- "[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?"

  safe_num <- function(x) {
    if (length(x) == 0 || is.na(x) || identical(x, "")) {
      return(NA_real_)
    }
    out <- suppressWarnings(as.numeric(x))
    ifelse(is.na(out), NA_real_, out)
  }

  extract_coefs_in_block <- function(block_lines) {
    if (length(block_lines) == 0) {
      return(tibble(term = character(), coef = numeric()))
    }

    model_idx <- grep("^\\s*Model:\\s*$", block_lines)
    if (length(model_idx) == 0) {
      return(tibble(term = character(), coef = numeric()))
    }

    coefs <- purrr::map_dfr(
      block_lines[(model_idx[1] + 1):length(block_lines)],
      function(line) {
        m <- stringr::str_match(
          line,
          paste0(
            "^\\s*([[:alnum:]_.]+)\\s+(",
            num_pat,
            ")\\s+(",
            num_pat,
            ")\\s+(",
            num_pat,
            ")(?:\\s+(",
            num_pat,
            "))?"
          )
        )
        if (all(is.na(m[1, ]))) {
          return(tibble(term = character(), coef = numeric()))
        }
        tibble(term = m[1, 2], coef = safe_num(m[1, 3]))
      }
    )
    coefs |>
      filter(!is.na(coef)) |>
      distinct(term, .keep_all = TRUE)
  }

  get_outcome_block <- function(lines, start_idx) {
    if (length(start_idx) == 0 || is.na(start_idx)) {
      return(character())
    }
    next_outcome <- grep(
      "^\\s*Outcome Variable:\\s*",
      lines,
      ignore.case = TRUE
    )
    next_outcome <- next_outcome[next_outcome > start_idx]
    end_idx <- if (length(next_outcome) > 0) {
      next_outcome[1] - 1
    } else {
      length(lines)
    }
    lines[start_idx:end_idx]
  }

  get_total_effect_block <- function(lines) {
    start <- grep(
      "^\\*+\\s*TOTAL EFFECT MODEL\\s*\\*+$",
      lines,
      ignore.case = TRUE
    )
    if (length(start) == 0) {
      return(character())
    }
    after <- lines[start[1]:length(lines)]
    end_rel <- grep("^\\*{10,}\\s*$", after)
    end_idx <- if (length(end_rel) > 1) {
      start[1] + end_rel[2] - 2
    } else {
      length(lines)
    }
    lines[start[1]:end_idx]
  }

  extract_effect_from_section <- function(lines, heading) {
    idx <- grep(heading, lines, ignore.case = TRUE)
    if (length(idx) == 0) {
      return(NA_real_)
    }
    look_ahead <- lines[(idx[1] + 1):min(length(lines), idx[1] + 8)]
    row <- grep(
      paste0("^\\s*(", num_pat, ")\\s+(", num_pat, ")\\s+(", num_pat, ")"),
      look_ahead,
      value = TRUE
    )
    if (length(row) == 0) {
      return(NA_real_)
    }
    m <- stringr::str_match(row[1], paste0("^\\s*(", num_pat, ")"))
    if (is.null(m) || nrow(m) == 0 || ncol(m) < 2) {
      return(NA_real_)
    }
    safe_num(m[1, 2])
  }

  lines <- stringr::str_split(output_text, "\\n")[[1]]

  # Parse n from model summary row using df1 + df2 + 1.
  n_df <- stringr::str_match(
    output_text,
    paste0(
      "(?m)^\\s*",
      num_pat,
      "\\s+",
      num_pat,
      "\\s+",
      num_pat,
      "\\s+",
      num_pat,
      "\\s+(",
      num_pat,
      ")\\s+(",
      num_pat,
      ")\\s+",
      num_pat,
      "\\s*$"
    )
  )
  n <- if (!is.null(n_df) && nrow(n_df) > 0) {
    safe_num(n_df[1, 2]) + safe_num(n_df[1, 3]) + 1
  } else {
    NA_real_
  }

  # Parse total/direct effect rows from the effects section.
  total_effect <- extract_effect_from_section(
    lines,
    "^\\s*Total effect of X on Y:\\s*$"
  )
  direct_effect <- extract_effect_from_section(
    lines,
    "^\\s*Direct effect of X on Y:\\s*$"
  )

  # Parse indirect table where BootULCI is printed on a separate line.
  ind_idx <- grep(
    "^\\s*Indirect effect\\(s\\) of X on Y:\\s*$",
    lines,
    ignore.case = TRUE
  )
  indirect_effect <- NA_real_
  indirect_llci <- NA_real_
  indirect_ulci <- NA_real_
  mediator_var <- NA_character_

  if (length(ind_idx) > 0) {
    ind_lines <- lines[(ind_idx[1] + 1):length(lines)]
    row1 <- grep(
      paste0(
        "^\\s*([[:alnum:]_.]+)\\s+(",
        num_pat,
        ")\\s+(",
        num_pat,
        ")\\s+(",
        num_pat,
        ")(?:\\s+(",
        num_pat,
        "))?\\s*$"
      ),
      ind_lines,
      value = TRUE
    )
    if (length(row1) > 0) {
      m1 <- stringr::str_match(
        row1[1],
        paste0(
          "^\\s*([[:alnum:]_.]+)\\s+(",
          num_pat,
          ")\\s+(",
          num_pat,
          ")\\s+(",
          num_pat,
          ")(?:\\s+(",
          num_pat,
          "))?\\s*$"
        )
      )
      mediator_var <- m1[1, 2]
      indirect_effect <- safe_num(m1[1, 3])
      indirect_llci <- safe_num(m1[1, 5])
      indirect_ulci <- safe_num(m1[1, 6])

      if (is.na(indirect_ulci)) {
        row2 <- grep(
          paste0("^\\s*", mediator_var, "\\s+(", num_pat, ")\\s*$"),
          ind_lines,
          value = TRUE
        )
        if (length(row2) > 0) {
          m2 <- stringr::str_match(
            row2[1],
            paste0("^\\s*", mediator_var, "\\s+(", num_pat, ")\\s*$")
          )
          indirect_ulci <- safe_num(m2[1, 2])
        }
      }
    }
  }

  # Parse a_path from first outcome block and b_path from second block.
  outcome_idx <- grep("^\\s*Outcome Variable:\\s*", lines, ignore.case = TRUE)
  med_block <- if (length(outcome_idx) >= 1) {
    get_outcome_block(lines, outcome_idx[1])
  } else {
    character()
  }
  y_block <- if (length(outcome_idx) >= 2) {
    get_outcome_block(lines, outcome_idx[2])
  } else {
    character()
  }

  med_coefs <- extract_coefs_in_block(med_block)
  y_coefs <- extract_coefs_in_block(y_block)

  total_block <- get_total_effect_block(lines)
  total_coefs <- extract_coefs_in_block(total_block)
  x_term <- total_coefs |>
    filter(term != "constant") |>
    slice(1) |>
    pull(term)

  if (length(x_term) > 0) {
    a_path <- med_coefs |>
      filter(term == x_term[1]) |>
      slice(1) |>
      pull(coef)
  } else {
    a_path <- numeric()
  }
  if (length(a_path) == 0) {
    a_path <- med_coefs |>
      filter(term != "constant") |>
      slice(1) |>
      pull(coef)
  }
  a_path <- if (length(a_path) == 0) NA_real_ else a_path

  b_path <- if (!is.na(mediator_var)) {
    y_coefs |>
      filter(term == mediator_var) |>
      slice(1) |>
      pull(coef)
  } else {
    numeric()
  }
  b_path <- if (length(b_path) == 0) NA_real_ else b_path

  tibble(
    n = n,
    total_effect = total_effect,
    direct_effect = direct_effect,
    indirect_effect = indirect_effect,
    indirect_llci = indirect_llci,
    indirect_ulci = indirect_ulci,
    a_path = a_path,
    b_path = b_path
  )
}

# --- tidy results ------------------------------------------------------

results <- models |>
  mutate(parsed = map(output, parse_process_output)) |>
  unnest(parsed) |>
  mutate(
    cov = map_chr(cov, ~ paste(.x, collapse = ", ")),
    icv_label = if_else(icv, "with icv", "no icv"),
    model_label = paste(dataset, structure, side, icv_label, sep = " | "),
    ci_excludes_zero = if_else(
      !is.na(indirect_llci) &
        !is.na(indirect_ulci) &
        ((indirect_llci > 0 & indirect_ulci > 0) |
          (indirect_llci < 0 & indirect_ulci < 0)),
      TRUE,
      FALSE,
      missing = FALSE
    )
  )

model_order <- results |>
  arrange(dataset, structure, side, icv_label) |>
  pull(model_label) |>
  unique()

results <- results |>
  mutate(model_label = factor(model_label, levels = model_order))

# --- forest plot -------------------------------------------------------

p_forest <- ggplot(
  results,
  aes(
    x = indirect_effect,
    y = model_label,
    xmin = indirect_llci,
    xmax = indirect_ulci
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0.2, na.rm = TRUE) +
  geom_point(aes(shape = ci_excludes_zero), size = 2.2, na.rm = TRUE) +
  facet_grid(structure ~ dataset, scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(
    x = "indirect effect (bootstrapped)",
    y = NULL,
    shape = "95% ci excludes 0"
  ) +
  theme_minimal(base_size = 11)

# --- specification grid ------------------------------------------------

spec_grid <- results |>
  select(model_label, dataset, structure, side, icv_label) |>
  pivot_longer(
    cols = c(dataset, structure, side, icv_label),
    names_to = "spec_type",
    values_to = "spec_value"
  ) |>
  mutate(
    model_label = factor(model_label, levels = levels(results$model_label)),
    spec_type = factor(
      spec_type,
      levels = c("dataset", "structure", "side", "icv_label")
    )
  )

p_specs <- ggplot(
  spec_grid,
  aes(x = spec_type, y = model_label, fill = spec_value)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  labs(x = NULL, y = NULL, fill = "spec") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_multiverse <- p_forest / p_specs + patchwork::plot_layout(heights = c(3, 2))

# --- save parsed results + figure --------------------------------------

outputs_dir <- file.path(proj_root, "outputs")
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  results,
  file.path(proj_root, "outputs", "process_multiverse_results.csv")
)

ggplot2::ggsave(
  filename = file.path(proj_root, "outputs", "process_multiverse_plot.png"),
  plot = p_multiverse,
  width = 14,
  height = 12,
  dpi = 300
)

# --- path diagrams -----------------------------------------------------

results_for_paths <- readr::read_csv(
  file.path(proj_root, "outputs", "process_multiverse_results.csv"),
  show_col_types = FALSE
)

node_template <- tibble(
  node = c("x", "m", "y"),
  node_label = c("x", "m", "y"),
  x = c(1, 2, 3),
  y = c(0, 1, 0)
)

nodes <- results_for_paths |>
  select(model_id, model_label) |>
  distinct() |>
  crossing(node_template)

edges <- bind_rows(
  results_for_paths |>
    transmute(
      model_id,
      model_label,
      path = "a",
      est = a_path,
      x = 1,
      y = 0,
      xend = 2,
      yend = 1
    ),
  results_for_paths |>
    transmute(
      model_id,
      model_label,
      path = "b",
      est = b_path,
      x = 2,
      y = 1,
      xend = 3,
      yend = 0
    ),
  results_for_paths |>
    transmute(
      model_id,
      model_label,
      path = "c'",
      est = direct_effect,
      x = 1,
      y = 0,
      xend = 3,
      yend = 0
    )
) |>
  mutate(
    edge_label = if_else(
      is.na(est),
      paste0(path, " = na"),
      paste0(path, " = ", sprintf("%.4f", est))
    ),
    label_x = (x + xend) / 2,
    label_y = (y + yend) / 2 + if_else(path == "c'", -0.18, 0.12)
  )

path_stats <- results_for_paths |>
  transmute(
    model_id,
    model_label,
    indirect_label = if_else(
      is.na(indirect_effect),
      "indirect = na",
      paste0(
        "indirect = ",
        sprintf("%.4f", indirect_effect),
        " [",
        sprintf("%.4f", indirect_llci),
        ", ",
        sprintf("%.4f", indirect_ulci),
        "]"
      )
    ),
    total_label = if_else(
      is.na(total_effect),
      "total = na",
      paste0("total = ", sprintf("%.4f", total_effect))
    ),
    ci_excludes_zero = ci_excludes_zero
  )

p_paths <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.09, "inches")),
    linewidth = 0.6,
    color = "grey20"
  ) +
  geom_text(
    data = edges,
    aes(x = label_x, y = label_y, label = edge_label),
    size = 2.8
  ) +
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = node_label),
    label.size = 0.2,
    fill = "white",
    size = 3
  ) +
  geom_text(
    data = path_stats,
    aes(
      x = 2,
      y = -0.45,
      label = total_label
    ),
    size = 2.8
  ) +
  geom_text(
    data = path_stats,
    aes(
      x = 2,
      y = -0.65,
      label = indirect_label,
      color = ci_excludes_zero
    ),
    size = 2.8
  ) +
  scale_color_manual(values = c(`TRUE` = "#0b7d20", `FALSE` = "#7a7a7a")) +
  facet_wrap(~model_label, ncol = 4) +
  coord_cartesian(xlim = c(0.7, 3.3), ylim = c(-0.8, 1.2), clip = "off") +
  theme_void(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    plot.margin = margin(8, 8, 8, 8)
  )

ggplot2::ggsave(
  filename = file.path(
    proj_root,
    "outputs",
    "process_multiverse_path_diagrams.png"
  ),
  plot = p_paths,
  width = 16,
  height = 14,
  dpi = 300
)

# --- html tabs ---------------------------------------------------------

make_tabs_html <- function(df, title = "process model outputs") {
  stopifnot(all(c("model_id", "output") %in% names(df)))

  tab_ids <- paste0("tab_", seq_len(nrow(df)))

  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }

  buttons <- Map(
    function(id, label, active) {
      tags$button(
        type = "button",
        class = paste("tab-btn", if (active) "active" else ""),
        `data-tab` = id,
        label
      )
    },
    tab_ids,
    df$model_id,
    seq_len(nrow(df)) == 1
  )

  panels <- Map(
    function(id, content, active) {
      tags$div(
        id = id,
        class = paste("tab-panel", if (active) "active" else ""),
        tags$pre(class = "code", HTML(esc(content)))
      )
    },
    tab_ids,
    df$output,
    seq_len(nrow(df)) == 1
  )

  tags$html(
    tags$head(
      tags$meta(charset = "utf-8"),
      tags$title(title),
      tags$style(HTML(
        "
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif; margin: 20px; }
        h1 { font-size: 18px; margin: 0 0 12px 0; }
        .tabs { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 12px; }
        .tab-btn {
          border: 1px solid #ccc; background: #f7f7f7; padding: 6px 10px; border-radius: 8px;
          cursor: pointer; font-size: 12px;
        }
        .tab-btn.active { background: #e9e9e9; border-color: #999; }
        .tab-panel { display: none; }
        .tab-panel.active { display: block; }
        pre.code {
          white-space: pre-wrap;
          background: #0b0b0b;
          color: #f2f2f2;
          padding: 12px;
          border-radius: 12px;
          border: 1px solid #222;
          overflow-x: auto;
          font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', monospace;
          font-size: 12px;
          line-height: 1.35;
        }
      "
      )),
      tags$script(HTML(
        "
        document.addEventListener('click', function(e) {
          if (!e.target.classList.contains('tab-btn')) return;
          const tabId = e.target.getAttribute('data-tab');
          document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
          document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
          e.target.classList.add('active');
          document.getElementById(tabId).classList.add('active');
        });
      "
      ))
    ),
    tags$body(
      tags$h1(title),
      tags$div(class = "tabs", buttons),
      tags$div(panels)
    )
  )
}

out_path <- file.path(proj_root, "outputs", "process_tabs.html")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

html <- make_tabs_html(
  models,
  title = "depressed-brain multiverse analysis: process outputs"
)
htmltools::save_html(html, file = out_path)

message("wrote: ", out_path)
utils::browseURL(out_path)

# --- sensitivity: nonharmonized + site as fixed effect ---------------

sens_spec <- spec |>
  filter(dataset == "nonharmonized") |>
  mutate(
    cov = map(
      icv,
      ~ if (.x) {
        c(cov_base, icv_var, "site_id_numeric")
      } else {
        c(cov_base, "site_id_numeric")
      }
    ),
    model_id = paste(
      "sensitivity",
      structure,
      side,
      ifelse(icv, "with_icv", "no_icv"),
      sep = "_"
    )
  )

nonharm_site <- nonharm_data |>
  mutate(site_id_numeric = as.numeric(factor(site)))

sens_models <- sens_spec |>
  mutate(
    output = pmap_chr(
      list(
        data = map(dataset, ~nonharm_site),
        cov = cov,
        y = rep(y_var, n()),
        x = rep(x_var, n()),
        m = m_var
      ),
      capture_process
    )
  )

sens_results <- sens_models |>
  mutate(parsed = map(output, parse_process_output)) |>
  unnest(parsed) |>
  mutate(
    cov = map_chr(cov, ~ paste(.x, collapse = ", ")),
    icv_label = if_else(icv, "with icv", "no icv"),
    model_label = paste("sensitivity", structure, side, icv_label, sep = " | "),
    ci_excludes_zero = if_else(
      !is.na(indirect_llci) &
        !is.na(indirect_ulci) &
        ((indirect_llci > 0 & indirect_ulci > 0) |
          (indirect_llci < 0 & indirect_ulci < 0)),
      TRUE,
      FALSE,
      missing = FALSE
    )
  )

sens_order <- sens_results |>
  arrange(structure, side, icv_label) |>
  pull(model_label) |>
  unique()

sens_results <- sens_results |>
  mutate(model_label = factor(model_label, levels = sens_order))

readr::write_csv(
  sens_results,
  file.path(proj_root, "outputs", "sensitivity_site_results.csv")
)

p_sens_forest <- ggplot(
  sens_results,
  aes(
    x = indirect_effect,
    y = model_label,
    xmin = indirect_llci,
    xmax = indirect_ulci
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0.2, na.rm = TRUE) +
  geom_point(aes(shape = ci_excludes_zero), size = 2.2, na.rm = TRUE) +
  facet_grid(structure ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(
    title = "sensitivity: nonharmonized + site_id as fixed effect",
    x = "indirect effect (bootstrapped)",
    y = NULL,
    shape = "95% ci excludes 0"
  ) +
  theme_minimal(base_size = 11)

ggplot2::ggsave(
  filename = file.path(proj_root, "outputs", "sensitivity_site_plot.png"),
  plot = p_sens_forest,
  width = 10,
  height = 8,
  dpi = 300
)

sens_html <- make_tabs_html(
  sens_models,
  title = "sensitivity analysis: nonharmonized + site_id fixed effect"
)
sens_out_path <- file.path(proj_root, "outputs", "sensitivity_site_tabs.html")
htmltools::save_html(sens_html, file = sens_out_path)

message("sensitivity results written: ", sens_out_path)
utils::browseURL(sens_out_path)
