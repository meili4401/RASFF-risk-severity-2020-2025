#!/usr/bin/env Rscript

# RASFF reviewer-response analysis
# Code 03: Figure 6A/B/C and Table 3 logistic model results

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(forcats)
  library(data.table)
  library(stringr)
  library(grid)
  library(forestploter)
})

analysis_dir <- "D:/桌面临时文件/已经发表论文/慧玲论文/RASFF食品大数据/重新分析/分析2"
data_file <- file.path(analysis_dir, "Final_RASFF_event_level_analysis_dataset.csv")

dt <- read_csv(data_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) %>%
  mutate(
    hazard_final = factor(
      hazard_final,
      levels = c(
        "Microbiological",
        "Labelling/Fraud/Quality",
        "Pesticides",
        "Mycotoxins",
        "Metals & pollutants",
        "Additives & composition",
        "Allergens",
        "Novel/GMO",
        "Physical/packaging",
        "Missing_hazard"
      )
    ),
    dist_group = factor(dist_group, levels = c("Multi-country", "Single-country", "Not specified"))
  )

extract_or <- function(model) {
  co <- summary(model)$coefficients
  res <- data.table(
    term = rownames(co),
    beta = co[, 1],
    se = co[, 2],
    p = co[, 4]
  )
  res <- res[term != "(Intercept)"]
  res[, aOR := exp(beta)]
  res[, CI_low := exp(beta - 1.96 * se)]
  res[, CI_high := exp(beta + 1.96 * se)]
  res[, p_fmt := fifelse(is.na(p), "", fifelse(p < 0.001, "<0.001", sprintf("%.4f", p)))]
  res[, OR_CI := sprintf("%.2f (%.2f-%.2f)", aOR, CI_low, CI_high)]
  res[]
}

prep_res_for_plot <- function(res) {
  res <- copy(res)
  res[, `Variable group` := fifelse(
    str_starts(term, "hazard_final"), "Hazard category",
    fifelse(str_starts(term, "dist_group"), "Distribution scope",
            fifelse(term == "year_c", "Year (per 1-year increase)", "Other"))
  )]
  res[, Level := term]
  res[str_starts(term, "hazard_final"), Level := str_replace(Level, "^hazard_final", "")]
  res[str_starts(term, "dist_group"), Level := str_replace(Level, "^dist_group", "")]
  res[term == "year_c", Level := "Year (per 1-year increase)"]
  res[, Level := str_squish(Level)]
  res[, `P value` := p_fmt]
  res[, `OR(95%CI)` := OR_CI]
  res[]
}

make_panel <- function(df, header_title, sort_by_or = TRUE) {
  df <- copy(df)
  if (sort_by_or) df <- df[order(-aOR)]
  header <- data.table(
    Level = header_title,
    ` ` = "",
    `OR(95%CI)` = "",
    `P value` = "",
    aOR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_
  )
  df[, ` ` := paste(rep(" ", 25), collapse = " ")]
  rbind(header, df[, .(Level, ` `, `OR(95%CI)`, `P value`, aOR, CI_low, CI_high)], fill = TRUE)
}

calc_device_height_px <- function(n_rows) max(1400, 420 + 120 * n_rows)
calc_pdf_height_in <- function(n_rows) max(5.8, 1.6 + 0.36 * n_rows)

tm <- forest_theme(
  base_size = 12,
  core = list(
    fg_params = list(hjust = 0, x = 0.02),
    bg_params = list(fill = c("#FFFFFF", "#F6F6F6"))
  ),
  colhead = list(
    fg_params = list(fontface = "bold", hjust = 0, x = 0.02)
  ),
  ci_pch = 16,
  ci_lwd = 2.2,
  ci_Theight = 0.2,
  refline_lwd = 1.2,
  refline_lty = "dashed",
  refline_col = "grey35"
)

plot_panel <- function(panel_df, outname, xlim, ticks_at,
                       xlab = "Adjusted OR (log scale)", x_trans = "log",
                       table_colwidths = c(0.60, 0.14, 0.18, 0.08)) {
  table_df <- as.data.frame(panel_df[, .(Level, ` `, `OR(95%CI)`, `P value`)])
  colnames(table_df) <- c("Level", " ", "OR(95%CI)", "P value")
  is_sum <- rep(FALSE, nrow(panel_df))
  is_sum[1] <- TRUE
  n_rows <- nrow(panel_df)

  png(file.path(analysis_dir, paste0(outname, ".png")),
      width = 5500, height = calc_device_height_px(n_rows), res = 600)
  grid.newpage()
  p <- forest(
    table_df,
    est = panel_df$aOR,
    lower = panel_df$CI_low,
    upper = panel_df$CI_high,
    ci_column = 2,
    is_summary = is_sum,
    ref_line = 1,
    xlim = xlim,
    ticks_at = ticks_at,
    xlab = xlab,
    x_trans = x_trans,
    colwidths = grid::unit(table_colwidths, "npc"),
    footnote = "",
    theme = tm
  )
  print(p)
  dev.off()

  pdf(file.path(analysis_dir, paste0(outname, ".pdf")), width = 9.5, height = calc_pdf_height_in(n_rows))
  grid.newpage()
  p2 <- forest(
    table_df,
    est = panel_df$aOR,
    lower = panel_df$CI_low,
    upper = panel_df$CI_high,
    ci_column = 2,
    is_summary = is_sum,
    ref_line = 1,
    xlim = xlim,
    ticks_at = ticks_at,
    xlab = xlab,
    x_trans = x_trans,
    colwidths = grid::unit(table_colwidths, "npc"),
    footnote = "",
    theme = tm
  )
  print(p2)
  dev.off()
}

# distribution scope may contribute to RASFF risk classification; therefore it is excluded from the primary model to avoid circular interpretation.
m_main <- glm(serious ~ year_c + hazard_final, family = binomial(), data = dt)
m_sens <- glm(serious ~ year_c + hazard_final + dist_group, family = binomial(), data = dt)
m_no_missing <- glm(
  serious ~ year_c + hazard_final,
  family = binomial(),
  data = dt %>% filter(hazard_final != "Missing_hazard") %>% mutate(hazard_final = fct_drop(hazard_final))
)

res_main_raw <- extract_or(m_main)
res_sens_raw <- extract_or(m_sens)
res_no_missing_raw <- extract_or(m_no_missing)

res_main <- prep_res_for_plot(res_main_raw)
res_sens <- prep_res_for_plot(res_sens_raw)

table3_main <- res_main_raw[, .(
  model = "Primary: no distribution scope",
  term, aOR, CI_low, CI_high, OR_CI, p, p_fmt
)]
table3_sens <- res_sens_raw[, .(
  model = "Sensitivity: with distribution scope",
  term, aOR, CI_low, CI_high, OR_CI, p, p_fmt
)]
table3_no_missing <- res_no_missing_raw[, .(
  model = "Sensitivity: excluding Missing_hazard",
  term, aOR, CI_low, CI_high, OR_CI, p, p_fmt
)]

table3 <- rbind(table3_main, table3_sens, table3_no_missing, fill = TRUE)
fwrite(table3, file.path(analysis_dir, "Table_3_Logistic_model_results.csv"))

df_main_hazard <- res_main[`Variable group` == "Hazard category",
                           .(Level, aOR, CI_low, CI_high, `OR(95%CI)`, `P value`)]
df_main_hazard_type <- df_main_hazard
fwrite(df_main_hazard_type, file.path(analysis_dir, "Table_3A_Primary_hazard_for_Figure6A.csv"))

panel_main <- make_panel(df_main_hazard_type, "Primary model: hazard category", TRUE)
plot_panel(
  panel_main,
  outname = "Figure_6A_Primary_Hazard_no_dist_group",
  xlim = c(0.05, 20),
  ticks_at = c(0.05, 0.1, 0.3, 1, 2, 5, 10, 20),
  xlab = "Adjusted OR (primary model; ref hazard = Microbiological)",
  table_colwidths = c(0.62, 0.12, 0.18, 0.08)
)

df_sens_dist <- res_sens[`Variable group` == "Distribution scope",
                         .(Level, aOR, CI_low, CI_high, `OR(95%CI)`, `P value`)]
fwrite(df_sens_dist, file.path(analysis_dir, "Table_3B_Sensitivity_distribution_scope_for_Figure6B.csv"))

panel_dist <- make_panel(df_sens_dist, "Sensitivity only: distribution scope (ref = Multi-country)", TRUE)
plot_panel(
  panel_dist,
  outname = "Figure_6B_Sensitivity_DistributionScope",
  xlim = c(0.2, 6),
  ticks_at = c(0.2, 0.3, 0.5, 1, 2, 3, 6),
  xlab = "Adjusted OR (sensitivity model only)",
  table_colwidths = c(0.62, 0.12, 0.18, 0.08)
)

year_compare <- rbind(
  data.table(model = "Primary: no distribution scope", res_main_raw[term == "year_c"]),
  data.table(model = "Sensitivity: with distribution scope", res_sens_raw[term == "year_c"]),
  data.table(model = "Sensitivity: excluding Missing_hazard", res_no_missing_raw[term == "year_c"]),
  fill = TRUE
)
year_panel_data <- year_compare[, .(
  Level = model,
  aOR,
  CI_low,
  CI_high,
  `OR(95%CI)` = OR_CI,
  `P value` = p_fmt
)]
fwrite(year_panel_data, file.path(analysis_dir, "Table_3C_Year_effect_comparison_for_Figure6C.csv"))

panel_year <- make_panel(year_panel_data, "Year trend: serious classification pattern", FALSE)
plot_panel(
  panel_year,
  outname = "Figure_6C_Year_Comparison",
  xlim = c(0.6, 1.2),
  ticks_at = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2),
  xlab = "Adjusted OR (per 1-year increase)",
  table_colwidths = c(0.62, 0.12, 0.18, 0.08)
)

cat("\nCode 03 complete: Figure 6A/B/C and Table 3.\n")
