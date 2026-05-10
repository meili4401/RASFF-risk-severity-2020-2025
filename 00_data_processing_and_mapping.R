#!/usr/bin/env Rscript

# RASFF reviewer-response analysis
# Code 00: data processing and mapping tables

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
})

analysis_dir <- "D:/桌面临时文件/已经发表论文/慧玲论文/RASFF食品大数据/重新分析/分析2"
input_file <- file.path(analysis_dir, "RASFF_window.csv")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

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

clean_chars <- function(x) {
  x <- str_squish(as.character(x))
  na_if(x, "")
}

first_official_hazard_category <- function(x) {
  x <- ifelse(is.na(x), "", x)
  matches <- str_extract_all(x, "\\{+[^{}]+\\}+")
  vapply(matches, function(m) {
    if (length(m) == 0) return(NA_character_)
    str_squish(str_remove_all(m[[1]], "[{}]"))
  }, character(1))
}

collapse_hazard <- function(category) {
  z <- str_to_lower(str_squish(category))
  case_when(
    is.na(z) | z == "" ~ "Missing_hazard",
    z %in% c(
      "pathogenic micro-organisms", "non-pathogenic micro-organisms",
      "biological contaminants", "parasitic infestation", "tses"
    ) ~ "Microbiological",
    z == "pesticide residues" ~ "Pesticides",
    z %in% c("mycotoxins", "natural toxins (other)") ~ "Mycotoxins",
    z %in% c(
      "heavy metals", "environmental pollutants", "industrial contaminants",
      "chemical contamination (other)",
      "residues of veterinary medicinal products", "radiation"
    ) ~ "Metals & pollutants",
    z %in% c(
      "food additives and flavourings", "composition", "feed additives"
    ) ~ "Additives & composition",
    z == "allergens" ~ "Allergens",
    z %in% c("novel food", "genetically modified", "gmo / novel food") ~ "Novel/GMO",
    z %in% c(
      "foreign bodies", "migration", "packaging defective / incorrect"
    ) ~ "Physical/packaging",
    z %in% c(
      "labelling absent/incomplete/incorrect", "adulteration / fraud",
      "organoleptic aspects", "not determined (other)",
      "poor or insufficient controls"
    ) ~ "Labelling/Fraud/Quality",
    TRUE ~ "Labelling/Fraud/Quality"
  )
}

make_dist_group <- function(distribution) {
  x <- clean_chars(distribution)
  n_countries <- ifelse(is.na(x), NA_integer_, lengths(str_split(x, "\\s*,\\s*")))
  case_when(
    is.na(x) ~ "Not specified",
    n_countries >= 2 ~ "Multi-country",
    n_countries == 1 ~ "Single-country",
    TRUE ~ "Not specified"
  )
}

raw <- read_csv(
  input_file,
  col_types = cols(.default = col_character()),
  locale = locale(encoding = "UTF-8"),
  name_repair = "minimal",
  progress = FALSE,
  show_col_types = FALSE
)
read_problems <- problems(raw)

raw_n <- nrow(raw)

cleaned <- raw %>%
  mutate(across(where(is.character), clean_chars)) %>%
  mutate(
    date_parsed = dmy_hms(date, quiet = TRUE),
    year = year(date_parsed),
    risk_decision = str_to_lower(risk_decision)
  )

window_pre_dedup <- cleaned %>%
  filter(year %in% 2020:2025, !is.na(reference))

food_window_pre_dedup <- window_pre_dedup %>%
  filter(str_to_lower(type) == "food")

analysis_data <- food_window_pre_dedup %>%
  arrange(reference, desc(date_parsed)) %>%
  distinct(reference, .keep_all = TRUE) %>%
  mutate(
    serious = if_else(risk_decision == "serious", 1L, 0L),
    class = if_else(serious == 1L, "Serious", "Non-serious"),
    hazard_official_first = first_official_hazard_category(hazards),
    hazard_final = collapse_hazard(hazard_official_first),
    hazard_final = factor(hazard_final, levels = hazard_levels),
    product_category = if_else(is.na(category), "Unknown_category", category),
    # distribution scope may contribute to RASFF risk classification; therefore it is excluded from the primary model to avoid circular interpretation.
    dist_group = make_dist_group(distribution),
    dist_group = factor(dist_group, levels = c("Multi-country", "Single-country", "Not specified")),
    year_c = year - 2020
  )

qc <- tibble(
  item = c(
    "raw_records",
    "readr_parse_problem_rows",
    "records_with_unparsed_date",
    "records_2020_2025_before_reference_deduplication",
    "non_food_records_excluded_before_analysis",
    "food_records_2020_2025_before_reference_deduplication",
    "duplicate_reference_records_removed",
    "event_level_records_after_reference_deduplication",
    "final_2020_2025_sample_size",
    "serious_n",
    "serious_prop"
  ),
  value = c(
    raw_n,
    nrow(read_problems),
    sum(is.na(cleaned$date_parsed)),
    nrow(window_pre_dedup),
    nrow(window_pre_dedup) - nrow(food_window_pre_dedup),
    nrow(food_window_pre_dedup),
    nrow(food_window_pre_dedup) - nrow(analysis_data),
    nrow(analysis_data),
    nrow(analysis_data),
    sum(analysis_data$serious),
    mean(analysis_data$serious)
  )
)

product_mapping <- analysis_data %>%
  distinct(original_category = product_category) %>%
  arrange(original_category) %>%
  mutate(
    mapped_product_group = original_category,
    mapping_note = "Original RASFF product category retained; no additional product recoding."
  )

hazard_mapping <- tibble(
  official_hazard_category = c(
    "pathogenic micro-organisms; non-pathogenic micro-organisms; biological contaminants; parasitic infestation; TSEs",
    "pesticide residues",
    "mycotoxins; natural toxins (other)",
    "heavy metals; environmental pollutants; industrial contaminants; chemical contamination (other); residues of veterinary medicinal products; radiation",
    "food additives and flavourings; composition; feed additives",
    "allergens",
    "novel food; genetically modified; GMO / novel food",
    "foreign bodies; migration; packaging defective / incorrect",
    "labelling absent/incomplete/incorrect; adulteration / fraud; organoleptic aspects; not determined (other); poor or insufficient controls",
    "missing hazards field"
  ),
  mapped_hazard_group = c(
    "Microbiological",
    "Pesticides",
    "Mycotoxins",
    "Metals & pollutants",
    "Additives & composition",
    "Allergens",
    "Novel/GMO",
    "Physical/packaging",
    "Labelling/Fraud/Quality",
    "Missing_hazard"
  ),
  mapping_note = c(
    rep("Collapsed from official RASFF hazard category extracted from braces in the hazards field.", 9),
    "Missing or blank hazards field retained as Missing_hazard."
  )
)

official_hazard_counts <- analysis_data %>%
  mutate(hazard_official_first = if_else(is.na(hazard_official_first), "Missing_hazard", hazard_official_first)) %>%
  count(hazard_official_first, sort = TRUE, name = "total")

processing_flow <- c(
  "# RASFF Data Processing Flow",
  "",
  "1. Read `RASFF_window.csv` from the `分析` folder using UTF-8 encoding and character columns.",
  "2. Trim leading/trailing whitespace in all character variables and convert empty strings to missing values.",
  "3. Parse `date` as day-month-year hour-minute-second and extract `year`.",
  "4. Keep records from 2020 to 2025 and remove records without `reference`.",
  "5. Deduplicate by `reference` to obtain event-level records.",
  "6. Define `serious = 1` only when `risk_decision == \"serious\"`; all other risk decisions are coded as 0.",
  "7. Extract the first official hazard category inside braces from `hazards` and collapse it to 10 mutually exclusive hazard groups.",
  "8. Retain blank or missing hazards as `Missing_hazard`.",
  "9. Retain original RASFF product `category` as `product_category`.",
  "10. Because this analysis focuses on food notifications, exclude records with `type != food` before event-level modelling.",
  "11. Define `dist_group` as `Multi-country`, `Single-country`, or `Not specified` for descriptive and sensitivity analyses only.",
  "12. Primary logistic model excludes distribution scope because it may contribute to RASFF serious classification and would risk circular interpretation.",
  "13. Interpret model results as RASFF serious classification patterns among food notifications, not true food safety risk."
)

write_csv(analysis_data, file.path(analysis_dir, "Final_RASFF_event_level_analysis_dataset.csv"))
write_csv(qc, file.path(analysis_dir, "QC_data_processing_summary.csv"))
write_csv(product_mapping, file.path(analysis_dir, "Appendix_Product_category_mapping.csv"))
write_csv(hazard_mapping, file.path(analysis_dir, "Appendix_Hazard_mapping_official_to_analysis_groups.csv"))
write_csv(official_hazard_counts, file.path(analysis_dir, "Appendix_Official_hazard_category_counts.csv"))
writeLines(processing_flow, file.path(analysis_dir, "Data_processing_flow.md"), useBytes = TRUE)

cat("\nCode 00 complete: data processing and mapping.\n")
cat("Final event-level sample size: ", nrow(analysis_data), "\n", sep = "")
