# missing data for 48 and 72 hours from  MecROX_Redox_Cleaned (1) which can be found in Mec_

#library
library(tidyverse)
library(openxlsx)
library(gtsummary) 

#phospholipid data 
pld_data = read.xlsx("data/mecrox/MecROX_Redox_Cleaned (1).xlsx", sheet = "All_Patients_All_TP", startRow =2)

    # clean phospholipid data
    pld_data = pld_data %>%
        mutate(Timepoint = as.numeric(str_remove_all(Timepoint, "[^0-9.]"))) %>% 
        mutate(across(-c(Patient_ID, Group), as.numeric))


# Phopsholipid percentage composition which was cleaned by myself  
add_pld =read.xlsx("data/mecrox/MecROX_pld_coc_percentage_cleaned.xlsx")

    #clean additional phospholipid data (keep only 48h and 72h timepoints)
    add_pld = add_pld %>%
        mutate(`Time.(h)` = as.numeric(str_remove_all(`Time.(h)`, "[^0-9.]"))) %>%
        filter(`Time.(h)` %in% c(48, 72)) %>%
        rename(Patient_ID = Patient.ID, Timepoint = `Time.(h)`)

# find shared numeric columns (excluding join keys)
shared_num <- intersect(
    names(pld_data)[sapply(pld_data, is.numeric)],
    names(add_pld)[sapply(add_pld, is.numeric)]
)
shared_num <- setdiff(shared_num, c("Patient_ID", "Timepoint"))

# rename shared columns in add_pld to avoid conflict
add_pld <- add_pld %>% rename_with(~ paste0(.x, "_add"), all_of(shared_num))

# save original patients
orig_patients <- pld_data %>% distinct(Patient_ID)

# join and coalesce
pld_data <- full_join(pld_data, add_pld, by = c("Patient_ID", "Timepoint"))
for (col in shared_num) {
    add_col <- paste0(col, "_add")
    if (add_col %in% names(pld_data)) {
        pld_data[[col]] <- coalesce(pld_data[[col]], pld_data[[add_col]])
        pld_data[[add_col]] <- NULL
    }
}

# remove patients not in original pld_data
pld_data <- pld_data %>% semi_join(orig_patients, by = "Patient_ID")

# remove if TFT_uM is NA and Timepoint is 48 or 72
pld_data <- pld_data %>% filter(!(is.na(TFT_uM) & Timepoint %in% c(48, 72)))            
