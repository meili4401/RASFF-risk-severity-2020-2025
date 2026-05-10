#!/usr/bin/env Rscript

# RASFF reviewer-response analysis
# Code 01: Figures 1-4 and Table 1 descriptive analysis

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

analysis_dir <- "D:/桌面临时文件/已经发表论文/慧玲论文/RASFF食品大数据/重新分析/分析2"
data_file <- file.path(analysis_dir, "Final_RASFF_event_level_analysis_dataset.csv")
nature_main <- "#3C5488FF"
nature_light <- "#D3D3D3"

dt <- read_csv(data_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))

base_theme <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold")
    )
}

desc_block <- function(data, group_var, section_name) {
  data %>%
    group_by({{ group_var }}) %>%
    summarise(
      total = n(),
      serious_n = sum(serious),
      serious_prop = serious_n / total,
      .groups = "drop"
    ) %>%
    rename(level = {{ group_var }}) %>%
    mutate(level = as.character(level)) %>%
    mutate(section = section_name, .before = 1) %>%
    arrange(section, desc(total))
}

table1 <- bind_rows(
  desc_block(dt, year, "Year"),
  desc_block(dt, hazard_final, "Hazard category"),
  desc_block(dt, dist_group, "Distribution scope"),
  desc_block(dt, product_category, "Product category")
) %>%
  mutate(serious_prop_percent = serious_prop * 100)

table1_top <- bind_rows(
  desc_block(dt, year, "Year"),
  desc_block(dt, hazard_final, "Hazard category"),
  desc_block(dt, dist_group, "Distribution scope"),
  desc_block(dt, product_category, "Top 15 product categories") %>% slice_head(n = 15)
) %>%
  mutate(serious_prop_percent = serious_prop * 100)

write_csv(table1, file.path(analysis_dir, "Table_1_Descriptive_analysis_full.csv"))
write_csv(table1_top, file.path(analysis_dir, "Table_1_Descriptive_analysis_main.csv"))

year_counts <- dt %>%
  count(year, name = "notifications") %>%
  arrange(year)

fig1 <- ggplot(year_counts, aes(x = factor(year), y = notifications)) +
  geom_col(fill = nature_main, color = "black", width = 0.75, linewidth = 0.4) +
  geom_text(aes(label = comma(notifications)), vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Year", y = "Number of RASFF notifications") +
  base_theme() +
  theme(panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5))

year_prop <- dt %>%
  mutate(class = factor(if_else(serious == 1L, "Serious", "Non-serious"), levels = c("Non-serious", "Serious"))) %>%
  count(year, class, name = "N") %>%
  group_by(year) %>%
  mutate(prop = N / sum(N), label = percent(prop, accuracy = 0.1)) %>%
  ungroup()

fig2 <- ggplot(year_prop, aes(x = factor(year), y = prop, fill = class)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4.3,
    fontface = "bold",
    color = "black"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = c("Serious" = nature_main, "Non-serious" = nature_light)) +
  labs(x = "Year", y = "Proportion of RASFF notifications", fill = NULL) +
  base_theme() +
  theme(
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 13)
  )

haz_all <- dt %>%
  count(hazard_final, name = "N") %>%
  arrange(N) %>%
  mutate(hazard_final = factor(hazard_final, levels = hazard_final))

fig3 <- ggplot(haz_all, aes(x = N, y = hazard_final)) +
  geom_col(fill = nature_main, color = "black", width = 0.7, linewidth = 0.25) +
  geom_text(aes(label = comma(N)), hjust = -0.1, size = 3.8, fontface = "bold") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of notifications", y = "Hazard category") +
  base_theme(15) +
  theme(panel.grid.major.x = element_line(color = "grey85", linewidth = 0.5))

cat_top <- dt %>%
  count(product_category, name = "N") %>%
  arrange(desc(N)) %>%
  slice_head(n = 15) %>%
  arrange(N) %>%
  mutate(product_category = factor(product_category, levels = product_category))

fig4 <- ggplot(cat_top, aes(x = N, y = product_category)) +
  geom_col(fill = nature_main, color = "black", width = 0.75, linewidth = 0.3) +
  geom_text(aes(label = comma(N)), hjust = -0.12, size = 4.0, fontface = "bold") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of notifications", y = NULL) +
  base_theme() +
  theme(panel.grid.major.x = element_line(color = "grey85", linewidth = 0.5))

ggsave(file.path(analysis_dir, "Figure_1_Annual_counts.png"), fig1, width = 8.5, height = 4.8, dpi = 600)
ggsave(file.path(analysis_dir, "Figure_1_Annual_counts.pdf"), fig1, width = 8.5, height = 4.8)
ggsave(file.path(analysis_dir, "Figure_2_Annual_serious_proportion.png"), fig2, width = 8.5, height = 4.8, dpi = 600)
ggsave(file.path(analysis_dir, "Figure_2_Annual_serious_proportion.pdf"), fig2, width = 8.5, height = 4.8)
ggsave(file.path(analysis_dir, "Figure_3_Hazard_distribution.png"), fig3, width = 9, height = 6, dpi = 600)
ggsave(file.path(analysis_dir, "Figure_3_Hazard_distribution.pdf"), fig3, width = 9, height = 6)
ggsave(file.path(analysis_dir, "Figure_4_Top15_product_categories.png"), fig4, width = 9, height = 7, dpi = 600)
ggsave(file.path(analysis_dir, "Figure_4_Top15_product_categories.pdf"), fig4, width = 9, height = 7)

cat("\nCode 01 complete: Figures 1-4 and Table 1.\n")
