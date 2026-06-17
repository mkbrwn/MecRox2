# R script to produce tables for the phospholipid data

source("src/0_clean_phopholipid.r")
library(huxtable)

dir.create("output/table/phospholipids", showWarnings = FALSE)

# identify numeric columns (excluding keys)
num_cols <- setdiff(
    names(pld_data)[sapply(pld_data, is.numeric)],
    c("Patient_ID", "Timepoint")
)

# produce a stratified table: one section per Timepoint, comparing Group
pld_table <- pld_data %>%
    mutate(Timepoint = paste0("h", Timepoint)) %>%
    tbl_strata(
        strata = Timepoint,
        ~ .x %>%
            select(Group, all_of(num_cols)) %>%
            tbl_summary(
                by = Group,
                statistic = all_continuous() ~ "{mean} ({sd})",
                digits = all_continuous() ~ 2,
                missing = "no"
            ) %>%
            add_difference(
                test = all_continuous() ~ "t.test",
                estimate_fun = all_continuous() ~ label_style_number(digits = 2)
            ) %>%
            add_n(col_label = "**N**") %>%
            modify_header(
                label ~ "**Phospholipid**",
                estimate ~ "**Mean Difference**"
            )
    )

# save outputs
prefix <- "output/table/phospholipids/phospholipid_by_group"
as_gt(pld_table) %>% gt::gtsave(paste0(prefix, ".png"))
as_hux_table(pld_table) %>% quick_xlsx(file = paste0(prefix, ".xlsx"))

message("Saved: ", prefix, ".png and .xlsx")
