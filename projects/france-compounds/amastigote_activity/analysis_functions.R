# analysis_functions.R  ────────────────────────────────────────────────
# One self-contained function: run_pipeline()
#   • cleans duplicated YOLO columns (e.g. "2" + "2.0")
#   • aggregates per replicate
#   • computes infection % and parasite-burden
#   • builds ggplots with Dunnett asterisks   (works with ≥1 dose + control)
#   • returns a list: rep_df, summary_df, plots (infection & burden)

run_pipeline <- function(csv_path = "amastigote_raw_c7c3c5c6.csv") {
  
  # ── libraries needed here ────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(readr);   library(dplyr);  library(tidyr);    library(stringr)
    library(ggplot2); library(scales); library(emmeans);  library(ggpubr)
  })
  
  # ── 1. read CSV ─────────────────────────────────────────────────────
  raw <- read_csv(csv_path, show_col_types = FALSE,
                  col_types = cols(.default = "i", image = "c"),
                  name_repair = "minimal")
  
  # ── 2. normalise duplicated YOLO columns (0 … 4, 0.0 … 4.0 …) ───────
  for (cls in 0:4) {
    pattern <- paste0("^", cls, "(\\.0+)?$")             # 2 2.0 2.00 …
    cols    <- grep(pattern, names(raw), value = TRUE)
    if (length(cols) == 0) next
    canon <- as.character(cls)
    if (!canon %in% names(raw)) raw[[canon]] <- 0L
    raw[[canon]] <- raw[[canon]] + rowSums(raw[cols], na.rm = TRUE)
    raw[setdiff(cols, canon)] <- NULL                    # drop extras
  }
  
  # ── 3. rename + split filename parts ────────────────────────────────
  df <- raw %>%
    rename(
      intracellular_amastigotes = `2`,
      infected_cells            = `1`,
      extracellular_amastigotes = `4`,
      uninfected_cells          = `3`,
      edge_cells                = `0`
    ) %>%
    separate(image, into = c("compound","conc","replicate","rest"),
             sep = "_", remove = FALSE, extra = "merge") %>%
    mutate(
      conc       = str_replace(conc, "(?i)um$",  " µM"),
      conc       = str_replace(conc, "(?i)^(cont|ctrl)$", "control"),
      replicate  = as.integer(replicate)
    )
  
  # helper: order factor levels control < low … high
  ordered_levels <- function(x) {
    ctrl  <- "control"
    other <- setdiff(unique(x), ctrl)
    nums  <- suppressWarnings(parse_number(other))
    other <- other[order(nums)]
    factor(x, levels = c(ctrl, other))
  }
  
  # ── 4. per-replicate aggregates ─────────────────────────────────────
  rep_df <- df %>%
    group_by(compound, conc, replicate) %>%
    summarise(across(
      c(intracellular_amastigotes, infected_cells,
        extracellular_amastigotes, uninfected_cells, edge_cells),
      sum, .names = "{.col}"), .groups = "drop") %>%
    mutate(
      total_macrophages = infected_cells + uninfected_cells,
      infection_pct     = if_else(total_macrophages > 0,
                                  infected_cells / total_macrophages, NA_real_),
      parasite_burden   = if_else(infected_cells > 0,
                                  intracellular_amastigotes / infected_cells,
                                  NA_real_),
      conc              = ordered_levels(conc)
    ) %>%
    filter(!is.na(conc) & conc != "")          # drop rows w/out concentration
  
  # ── 5. summary table (mean ± sd) ────────────────────────────────────
  summary_df <- rep_df %>%
    group_by(compound, conc) %>%
    summarise(
      infection_mean = mean(infection_pct,   na.rm = TRUE),
      infection_sd   = sd  (infection_pct,   na.rm = TRUE),
      burden_mean    = mean(parasite_burden, na.rm = TRUE),
      burden_sd      = sd  (parasite_burden, na.rm = TRUE),
      .groups = "drop")
  
  # ── 6. helper for Dunnett labels ────────────────────────────────────
  get_dunnett_labels <- function(dat, yvar) {
    dat_ok <- dat[!is.na(dat[[yvar]]), ]
    if (!any(as.character(dat_ok$conc) == "control"))      return(data.frame())
    if (nlevels(droplevels(dat_ok$conc)) < 2)              return(data.frame())
    
    mod <- aov(reformulate("conc", yvar), data = dat_ok)
    cmp <- contrast(emmeans(mod, "conc"),
                    "trt.vs.ctrl", ref = "control", adjust = "dunnett") |>
      summary()
    
    p2symb <- function(p) {
      if (is.na(p) | p >= .05) "" else
        if (p < .0001) "****" else
          if (p < .001) "***" else
            if (p < .01) "**" else "*"
    }
    
    data.frame(
      group1   = "control",
      group2   = sub(" - .*", "", cmp$contrast),   # dose only
      p.signif = vapply(cmp$p.value, p2symb, character(1)),
      y.position =
        max(dat_ok[[yvar]], na.rm = TRUE) * 1.05 +
        seq(0,
            by = diff(range(dat_ok[[yvar]], na.rm = TRUE))*0.05,
            length.out = nrow(cmp))
    ) |>
      filter(p.signif != "")
  }
  
  # ── 7. build plots (skip metric if all NA) ──────────────────────────
  plots <- list()
  for (metric in c("infection_pct","parasite_burden")) {
    plot_df <- rep_df %>% filter(!is.na(.data[[metric]]))
    if (nrow(plot_df) == 0) { plots[[metric]] <- NULL; next }
    
    ylab  <- if (metric == "infection_pct")
      "Infection (%) of macrophages"
    else
      "Parasite burden (amastigotes / infected cell)"
    
    scale_y <- if (metric == "infection_pct")
      scale_y_continuous(labels = percent_format(accuracy = 1))
    else
      scale_y_continuous()
    
    p <- ggplot(plot_df, aes(conc, .data[[metric]], fill = conc)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
      facet_wrap(~ compound, scales = "free_x") +
      scale_y +
      labs(x = "Concentration", y = ylab) +
      theme_bw() +
      theme(legend.position = "none")
    
    lab_df <- plot_df %>%
      group_by(compound) %>%
      group_modify(~ get_dunnett_labels(.x, metric))
    
    if (nrow(lab_df) > 0)
      p <- p + ggpubr::stat_pvalue_manual(
        lab_df, label = "p.signif",
        hide.ns = TRUE,
        tip.length = 0.01,
        step.increase = 0.02,
        inherit.aes  = FALSE)
    
    plots[[metric]] <- p
  }
  
  # ── 8. return list ---------------------------------------------------
  list(rep_df = rep_df,
       summary_df = summary_df,
       plots = plots)
}
