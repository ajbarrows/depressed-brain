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

pacman::p_load("arrow", "dplyr", "tidyr", "purrr", "tibble", "htmltools")

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
