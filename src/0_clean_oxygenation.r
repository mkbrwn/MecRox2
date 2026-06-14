# script to clean the meta data for the oxygenation data

#library 
library(tidyverse)
library(openxlsx)
library(gtsummary) 

#load data 
oxygenation_data <- read.xlsx("data/meta_oxygenation/UKRoxData_V5.xlsx")

#select only the columns that are needed for the analysis
o2_data = oxygenation_data %>%
    select( MECROXStudy, IMVStart, UKRoxTime, SpO2Time, SpO2Value, Treatment,
            `PFRatio_Between_IMVStartTime_&_UKRoxTime`, `PFRatio_Between_UKRoxTime_+_5DaysUKRoxTime`,
            `FiO2_IMVStart_to_UKRox_Date`, `FiO2_UKRox_to_UKRox+5Days_Date`,
            `PO2_IMVStart_to_UKRox_Date`, `PO2_UKRox_to_UKRox+5Days_Date`, `UKRox_to_UKRox+5Days_Date./.time.of.last.Intubation/extubation?`
            ) %>% 
    filter( !is.na(MECROXStudy), !is.na(SpO2Time), !is.na(SpO2Value)) %>% 
    mutate( extubation = `UKRox_to_UKRox+5Days_Date./.time.of.last.Intubation/extubation?`) %>%
    mutate( SpO2Value = as.numeric(SpO2Value),
            SpO2Time  = convertToDateTime(SpO2Time))

# combine paired columns (paste both time periods) and remove originals
o2_data <- o2_data %>%
    mutate(
        PFRatio    = coalesce(str_c(`PFRatio_Between_IMVStartTime_&_UKRoxTime`, `PFRatio_Between_UKRoxTime_+_5DaysUKRoxTime`, sep = ", "), `PFRatio_Between_IMVStartTime_&_UKRoxTime`, `PFRatio_Between_UKRoxTime_+_5DaysUKRoxTime`),
        FiO2       = coalesce(str_c(`FiO2_IMVStart_to_UKRox_Date`, `FiO2_UKRox_to_UKRox+5Days_Date`, sep = ", "), `FiO2_IMVStart_to_UKRox_Date`, `FiO2_UKRox_to_UKRox+5Days_Date`),
        PO2        = coalesce(str_c(`PO2_IMVStart_to_UKRox_Date`, `PO2_UKRox_to_UKRox+5Days_Date`, sep = ", "), `PO2_IMVStart_to_UKRox_Date`, `PO2_UKRox_to_UKRox+5Days_Date`)
    ) %>%
    select(-`PFRatio_Between_IMVStartTime_&_UKRoxTime`, -`PFRatio_Between_UKRoxTime_+_5DaysUKRoxTime`,
           -`FiO2_IMVStart_to_UKRox_Date`, -`FiO2_UKRox_to_UKRox+5Days_Date`,
           -`PO2_IMVStart_to_UKRox_Date`, -`PO2_UKRox_to_UKRox+5Days_Date`,-`UKRox_to_UKRox+5Days_Date./.time.of.last.Intubation/extubation?`)

# parse "(dd/mm/yyyy HH:MM)value" strings into long format
# helper: splits a cell into individual entries, extracts datetime and value
parse_time_value <- function(x) {
    entries <- str_extract_all(x, "\\(\\d{2}/\\d{2}/\\d{4} \\d{2}:\\d{2}\\)[0-9.]+")[[1]]
    datetime <- str_extract(entries, "(?<=\\()[^\\)]+")
    value    <- as.numeric(str_extract(entries, "(?<=\\))[0-9.]+"))
    tibble(datetime = datetime, value = value)
}

o2_data_long = o2_data %>%
    mutate(UKRoxTime = ymd_hm(UKRoxTime)) %>%
    select(MECROXStudy, IMVStart, UKRoxTime, Treatment, SpO2Time, SpO2Value, PFRatio, FiO2, PO2) %>%
    pivot_longer(cols = c(PFRatio, FiO2, PO2), names_to = "variable", values_to = "raw") %>%
    filter(!is.na(raw), raw != "") %>%
    mutate(parsed = map(raw, parse_time_value)) %>%
    unnest(parsed) %>%
    mutate(datetime = dmy_hm(datetime)) %>%
    select(-raw)

# add SpO2 (already in wide format: one datetime + value per patient)
spo2_long <- o2_data %>%
    mutate(UKRoxTime = ymd_hm(UKRoxTime)) %>%
    select(MECROXStudy, IMVStart, UKRoxTime, Treatment, SpO2Time, SpO2Value) %>%
    mutate(variable = "SpO2", datetime = SpO2Time, value = SpO2Value) %>%
    select(-SpO2Time, -SpO2Value) %>%
    filter(!is.na(value))

o2_data_long <- bind_rows(o2_data_long, spo2_long)

o2_data_long = o2_data_long %>%
    arrange(MECROXStudy, variable, datetime) %>%
    select( -SpO2Time, -SpO2Value) %>%
    filter(!(variable == "SpO2" & value < 80)) 

# calculate variables required for analysis 
o2_data_long = o2_data_long %>%
    mutate(Treatment = recode(Treatment, "Usual" = "Usual care")) %>%
    mutate(TimeSinceRandomisation = as.numeric(difftime(datetime, UKRoxTime, units = "hours"))) %>%
    group_by(MECROXStudy) %>% mutate(Maxtime = max(TimeSinceRandomisation, na.rm = T))

# filter if time since randomisation is greater than 120 and < -12 (eligibility criteria is 12 hours pre-randomisation)
o2_data_long = o2_data_long %>%
    filter(TimeSinceRandomisation <= 121 & TimeSinceRandomisation >= -12)

# deduplicate rows (if any) based on MECROXStudy, variable, and datetime
o2_data_long = o2_data_long %>%
    distinct(MECROXStudy, variable, datetime, .keep_all = TRUE)

#filter rows 



# save data as rds file
saveRDS(o2_data_long, "data/processed/cleaned_oxygenation_data.rds")
# saveRDS(o2_data_long, "data/processed/cleaned_oxygenation_data.rds")
