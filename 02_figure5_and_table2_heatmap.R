#!/usr/bin/env Rscript

# RASFF reviewer-response analysis
# Code 02: Figure 5 hazard x top 15 product category heatmap and Table 2

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

analysis_dir <- "D:/桌面临时文件/已经发表论文/慧玲论文/RASFF食品大数据/重新分析/分析2"
data_file <- file.path(analysis_dir, "Final_RASFF_event_level_analysis_dataset.csv")
nature_main <- "#3C5488FF"

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

dt <- read_csv(data_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))

top_cats <- dt %>%
  count(product_category, sort = TRUE) %>%
  slice_head(n = 15) %>%
  pull(product_category)

table2_long <- dt %>%
  filter(product_category %in% top_cats) %>%
  count(hazard_final, product_category, name = "total") %>%
  arrange(hazard_final, product_category)

table2_fullgrid <- expand_grid(
  hazard_final = hazard_levels,
  product_category = top_cats
) %>%
  left_join(table2_long, by = c("hazard_final", "product_category")) %>%
  mutate(
    total = replace_na(total, 0L),
    log10_count_plus_1 = log10(total + 1)
  )

table2_matrix <- table2_fullgrid %>%
  select(hazard_final, product_category, total) %>%
  pivot_wider(names_from = product_category, values_from = total)

write_csv(table2_long, file.path(analysis_dir, "Table_2_Hazard_by_top15_product_counts_long.csv"))
write_csv(table2_fullgrid, file.path(analysis_dir, "Table_2_Hazard_by_top15_product_counts_fullgrid.csv"))
write_csv(table2_matrix, file.path(analysis_dir, "Table_2_Hazard_by_top15_product_counts_matrix.csv"))

plot_data <- table2_fullgrid %>%
  mutate(
    hazard_final = factor(hazard_final, levels = rev(hazard_levels)),
    product_category = factor(product_category, levels = top_cats)
  )

fig5 <- ggplot(plot_data, aes(x = product_category, y = hazard_final, fill = log10_count_plus_1)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = nature_main, name = "log10(Count+1)") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(analysis_dir, "Figure_5_Hazard_top15_product_heatmap.png"), fig5, width = 14, height = 7, dpi = 600)
ggsave(file.path(analysis_dir, "Figure_5_Hazard_top15_product_heatmap.pdf"), fig5, width = 14, height = 7)

cat("\nCode 02 complete: Figure 5 and Table 2.\n")
