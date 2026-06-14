# produce summary tables for oxygenation metrics

library(tidyverse)
library(gtsummary)
library(huxtable)

# load cleaned data (run cleaning script if needed)
if (!file.exists("data/processed/cleaned_oxygenation_data.rds")) {
    source("src/0_clean_oxygenation.r")
}
o2_data_long <- readRDS("data/processed/cleaned_oxygenation_data.rds")

dir.create("output/table/oxygenation_metrics", recursive = TRUE, showWarnings = FALSE)

# --- helper: produce patient-level and observation-level tables ---------------
make_tables <- function(data, var_name) {

    dat <- data %>% filter(variable == var_name)

    # --- patient-level: mean per patient in each 12h window --------------------
    patient_summary <- dat %>%
        filter(!is.na(TimeSinceRandomisation), !is.na(value)) %>%
        mutate(TimeWindow = ifelse(
            TimeSinceRandomisation < 0, -12,
            floor(TimeSinceRandomisation / 12) * 12
        )) %>%
        filter(TimeWindow <= 72) %>%
        group_by(MECROXStudy, TimeWindow) %>%
        summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(
            id_cols    = MECROXStudy,
            names_from = TimeWindow,
            values_from = mean_val,
            names_glue = "{TimeWindow}-{TimeWindow + 12}"
        ) %>%
        left_join(dat %>% distinct(MECROXStudy, Treatment), by = "MECROXStudy") %>%
        select(-MECROXStudy)

    tbl_patient <- patient_summary %>%
        tbl_summary(
            by         = Treatment,
            statistic  = list(all_continuous() ~ "{mean} ({sd})"),
            digits     = all_continuous() ~ 2,
            missing    = "no"
        ) %>%
        add_difference(
            test = all_continuous() ~ "t.test",
            estimate_fun = all_continuous() ~ label_style_number(digits = 2)
        ) %>%
        add_n(col_label = "**Patients**") %>%
        modify_header(
            label    ~ "**Time Since Randomisation (hours)**",
            estimate ~ "**Mean Difference**"
        )

    # save patient-level
    prefix <- paste0("output/table/oxygenation_metrics/", tolower(var_name))
    as_gt(tbl_patient) %>% gt::gtsave(paste0(prefix, "_patient.png"))
    as_hux_table(tbl_patient) %>% quick_xlsx(file = paste0(prefix, "_patient.xlsx"))

    # --- observation-level: each individual reading in each 12h window ---------
    obs_summary <- dat %>%
        filter(!is.na(TimeSinceRandomisation), !is.na(value)) %>%
        mutate(TimeWindow = ifelse(
            TimeSinceRandomisation < 0, -12,
            floor(TimeSinceRandomisation / 12) * 12
        )) %>%
        filter(TimeWindow <= 72) %>%
        ungroup() %>%
        mutate(obs_id = row_number()) %>%
        pivot_wider(
            id_cols     = c(obs_id, Treatment),
            names_from  = TimeWindow,
            values_from = value,
            names_glue = "{TimeWindow}-{TimeWindow + 12}"
        ) %>%
        select(-obs_id) %>%
        rename_with(~ gsub("^h", "", .x), starts_with("h"))

    tbl_obs <- obs_summary %>%
        tbl_summary(
            by         = Treatment,
            statistic  = list(all_continuous() ~ "{mean} ({sd})"),
            digits     = all_continuous() ~ 2,
            missing    = "no"
        ) %>%
        add_difference(
            test = all_continuous() ~ "t.test",
            estimate_fun = all_continuous() ~ label_style_number(digits = 2)
        ) %>%
        add_n(col_label = "**Observations**") %>%
        modify_header(
            label    ~ "**Time Since Randomisation (hours)**",
            estimate ~ "**Mean Difference**"
        )

    # save observation-level
    as_gt(tbl_obs) %>% gt::gtsave(paste0(prefix, "_obs.png"))
    as_hux_table(tbl_obs) %>% quick_xlsx(file = paste0(prefix, "_obs.xlsx"))

    message("Saved tables for: ", var_name)
}

# --- generate tables for each variable ----------------------------------------

for (v in c("PFRatio", "SpO2", "PO2", "FiO2")) {
    make_tables(o2_data_long, v)
}

message("All tables saved to output/table/oxygenation_metrics/")