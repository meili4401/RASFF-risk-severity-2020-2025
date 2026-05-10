#!/usr/bin/env Rscript

# Missing_hazard sensitivity analysis aligned with the primary model.
# Purpose: assess whether excluding Missing_hazard changes the primary hazard associations.

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
out_dir <- file.path(analysis_dir, "敏感分析去除缺失")
data_file <- file.path(analysis_dir, "Final_RASFF_event_level_analysis_dataset.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hazard_levels <- c(
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

prep_hazard_res <- function(res, model_name) {
  copy(res)[str_starts(term, "hazard_final")] %>%
    as_tibble() %>%
    mutate(
      model = model_name,
      hazard = str_squish(str_replace(term, "^hazard_final", ""))
    ) %>%
    select(model, hazard, term, aOR, CI_low, CI_high, OR_CI, p, p_fmt)
}

make_panel <- function(df, header_title) {
  header <- data.table(
    Level = header_title,
    ` ` = "",
    `OR(95%CI)` = "",
    `P value` = "",
    aOR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_
  )
  df <- copy(df)
  df[, ` ` := paste(rep(" ", 25), collapse = " ")]
  rbind(header, df[, .(Level, ` `, `OR(95%CI)`, `P value`, aOR, CI_low, CI_high)], fill = TRUE)
}

calc_device_height_px <- function(n_rows) max(1500, 420 + 120 * n_rows)
calc_pdf_height_in <- function(n_rows) max(6, 1.6 + 0.36 * n_rows)

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
                       xlab = "Adjusted OR (log scale)",
                       table_colwidths = c(0.62, 0.12, 0.18, 0.08)) {
  table_df <- as.data.frame(panel_df[, .(Level, ` `, `OR(95%CI)`, `P value`)])
  colnames(table_df) <- c("Level", " ", "OR(95%CI)", "P value")
  is_sum <- rep(FALSE, nrow(panel_df))
  is_sum[1] <- TRUE

  png(file.path(out_dir, paste0(outname, ".png")),
      width = 5500, height = calc_device_height_px(nrow(panel_df)), res = 600)
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
    x_trans = "log",
    colwidths = grid::unit(table_colwidths, "npc"),
    footnote = "",
    theme = tm
  )
  print(p)
  dev.off()

  pdf(file.path(out_dir, paste0(outname, ".pdf")), width = 9.5, height = calc_pdf_height_in(nrow(panel_df)))
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
    x_trans = "log",
    colwidths = grid::unit(table_colwidths, "npc"),
    footnote = "",
    theme = tm
  )
  print(p2)
  dev.off()
}

dt <- read_csv(data_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8")) %>%
  mutate(
    hazard_final = factor(hazard_final, levels = hazard_levels),
    dist_group = factor(dist_group, levels = c("Multi-country", "Single-country", "Not specified"))
  )

dt_no_missing_hazard <- dt %>%
  filter(hazard_final != "Missing_hazard") %>%
  mutate(hazard_final = fct_drop(hazard_final))

# Primary model from the main analysis.
m_primary <- glm(serious ~ year_c + hazard_final, family = binomial(), data = dt)

# Missing-value sensitivity model: same formula, excluding Missing_hazard records.
m_excluding_missing <- glm(
  serious ~ year_c + hazard_final,
  family = binomial(),
  data = dt_no_missing_hazard
)

primary_res <- extract_or(m_primary)
excluding_missing_res <- extract_or(m_excluding_missing)

primary_hazard <- prep_hazard_res(primary_res, "Primary: Missing_hazard retained")
excluding_missing_hazard <- prep_hazard_res(excluding_missing_res, "Sensitivity: Missing_hazard excluded")

hazard_comparison <- full_join(
  primary_hazard %>%
    select(hazard, primary_OR = aOR, primary_CI_low = CI_low, primary_CI_high = CI_high,
           primary_OR_95CI = OR_CI, primary_p = p, primary_p_fmt = p_fmt),
  excluding_missing_hazard %>%
    select(hazard, sensitivity_OR = aOR, sensitivity_CI_low = CI_low, sensitivity_CI_high = CI_high,
           sensitivity_OR_95CI = OR_CI, sensitivity_p = p, sensitivity_p_fmt = p_fmt),
  by = "hazard"
) %>%
  mutate(
    OR_ratio_sensitivity_vs_primary = sensitivity_OR / primary_OR,
    interpretation = case_when(
      is.na(sensitivity_OR) ~ "Not estimable after excluding Missing_hazard",
      abs(log(OR_ratio_sensitivity_vs_primary)) < log(1.10) ~ "Very similar",
      abs(log(OR_ratio_sensitivity_vs_primary)) < log(1.25) ~ "Similar",
      TRUE ~ "Meaningfully changed"
    )
  ) %>%
  arrange(match(hazard, hazard_levels))

write_csv(
  tibble(
    item = c(
      "food_only_event_level_records",
      "records_with_Missing_hazard",
      "records_after_excluding_Missing_hazard",
      "primary_formula",
      "missing_value_sensitivity_formula",
      "main_comparison"
    ),
    value = c(
      nrow(dt),
      sum(dt$hazard_final == "Missing_hazard", na.rm = TRUE),
      nrow(dt_no_missing_hazard),
      "serious ~ year_c + hazard_final",
      "serious ~ year_c + hazard_final, excluding Missing_hazard",
      "Compare hazard-specific ORs before and after excluding Missing_hazard"
    )
  ),
  file.path(out_dir, "QC_missing_hazard_OR_comparison.csv")
)

write_csv(primary_hazard, file.path(out_dir, "Table_MH1_Primary_hazard_ORs.csv"))
write_csv(excluding_missing_hazard, file.path(out_dir, "Table_MH2_Excluding_Missing_hazard_ORs.csv"))
write_csv(hazard_comparison, file.path(out_dir, "Table_MH3_Hazard_OR_comparison_primary_vs_excluding_missing.csv"))

plot_rows <- bind_rows(primary_hazard, excluding_missing_hazard) %>%
  filter(hazard != "Missing_hazard") %>%
  mutate(
    hazard = factor(hazard, levels = rev(hazard_levels[hazard_levels != "Missing_hazard"])),
    model = factor(model, levels = c("Primary: Missing_hazard retained", "Sensitivity: Missing_hazard excluded")),
    Level = paste0(as.character(hazard), " - ", as.character(model)),
    `OR(95%CI)` = OR_CI,
    `P value` = p_fmt
  ) %>%
  arrange(hazard, model) %>%
  as.data.table()

panel <- make_panel(
  plot_rows[, .(Level, aOR, CI_low, CI_high, `OR(95%CI)`, `P value`)],
  "Missing_hazard sensitivity: hazard OR comparison"
)

plot_panel(
  panel,
  outname = "Figure_MH1_Hazard_OR_comparison_primary_vs_excluding_missing",
  xlim = c(0.05, 20),
  ticks_at = c(0.05, 0.1, 0.3, 1, 2, 5, 10, 20),
  xlab = "Adjusted OR (ref hazard = Microbiological)"
)

cat("\nMissing_hazard OR comparison sensitivity analysis complete.\n")
cat("Output folder: ", out_dir, "\n", sep = "")
cat("Primary records: ", nrow(dt), "\n", sep = "")
cat("Missing_hazard records excluded in sensitivity model: ", sum(dt$hazard_final == "Missing_hazard", na.rm = TRUE), "\n", sep = "")
cat("Sensitivity records: ", nrow(dt_no_missing_hazard), "\n", sep = "")
