# R script to produce figures for phospholipid data

source("src/0_clean_phopholipid.r")


# --- Figures: change from baseline boxplots -----------------------------------

dir.create("output/figures", showWarnings = FALSE)

plot_change_from_baseline <- function(data, cols, panel_title, filename) {

    # calculate change from baseline (timepoint 0) for each patient
    baseline <- data %>%
        filter(Timepoint == 0) %>%
        select(Patient_ID, Group, all_of(cols)) %>%
        pivot_longer(cols = all_of(cols), names_to = "variable", values_to = "baseline")

    change_data <- data %>%
        filter(Timepoint != 0) %>%
        select(Patient_ID, Group, Timepoint, all_of(cols)) %>%
        pivot_longer(cols = all_of(cols), names_to = "variable", values_to = "value") %>%
        left_join(baseline, by = c("Patient_ID", "Group", "variable")) %>%
        mutate(change = value - baseline,
               Timepoint = factor(paste0("h", Timepoint)))

    # compute medians per Timepoint × variable × Group for connecting line
    medians <- change_data %>%
        group_by(Timepoint, variable, Group) %>%
        summarise(median_change = median(change, na.rm = TRUE), .groups = "drop")

    # compute p-values: paired wilcoxon test of change vs 0 at each timepoint
    pvals <- change_data %>%
        group_by(variable, Timepoint, Group) %>%
        summarise(
            p = wilcox.test(change, mu = 0)$p.value,
            med = median(change, na.rm = TRUE),
            y_max = max(change, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(
            label = ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**",
                     ifelse(p < 0.05, "*", "ns"))),
            x_pos = case_when(
                Group == "Conservative" & Timepoint == "h48" ~ 0.78,
                Group == "Usual"        & Timepoint == "h48" ~ 1.22,
                Group == "Conservative" & Timepoint == "h72" ~ 1.78,
                TRUE                                         ~ 2.22
            ),
            # place p-value just above the max value for that group/timepoint
            y_pos = y_max * 1.1
        )

    # labels
    var_labels <- setNames(
        str_replace_all(cols, "_", " "),
        cols
    )

    p <- ggplot(change_data, aes(x = Timepoint, y = change, fill = Group)) +
        geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, size = 1, alpha = 0.4, aes(colour = Group)) +
        geom_line(data = medians, aes(x = Timepoint, y = median_change,
                                       group = Group, colour = Group),
                  linewidth = 1) +
        geom_point(data = medians, aes(x = Timepoint, y = median_change,
                                        colour = Group), size = 3) +
        geom_text(data = pvals, aes(x = x_pos, y = y_pos, label = label),
                  size = 4, fontface = "bold", colour = "grey20") +
        facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = var_labels)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        labs(x = "Timepoint", y = "Change from baseline", title = panel_title,
             subtitle = "Wilcoxon signed-rank test: change vs 0") +
        scale_colour_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        theme_bw(base_size = 14) +
        theme(
            strip.text = element_text(face = "bold", size = 11),
            strip.background = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 10, colour = "grey40"),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12)
        )

    ggsave(paste0("output/figures/", filename), plot = p, width = 12, height = 8, dpi = 300)
    message("Saved: ", filename)
}

# Specifying the columns for oxidative stress and inorganic species, free thiols, total thiols, and phospholipids
oxidative_stress_cols <- c("TFT_uM", "TFT_umol_g_protein", "Protein_g_L",
                            "FRAP_uM", "TBARS_nM", "RXNO_nM",
                            "Nitrite_uM", "Nitrate_uM", "Thiosulfate_nM", "Sulfate_uM")

Free_thiols_nM <- c( "GSH_free_nM","GSSG_free_nM","CyS_free_nM", "CySS_free_nM","HCyS_free_nM", "HCySS_free_nM", "Cysteamine_free_nM","GluCyS_free_nM", "CysGly_free_nM", "NAC_free_nM", "Sulfide_free_nM"
)

Total_thiols_nM <- c("GSH_total_uM", "CyS_total_uM", "HCyS_total_uM", "Cysteamine_total_uM", "GluCyS_total_uM", "CysGly_total_uM", "NAC_total_uM", "Sulfide_total_uM" 
)

Phospholipids_perc_competition_PC <- c( 

  "PC16:0/14:0", "PC16:0a/16:0", "PC16:0/16:1", "PC16:0/16:0",
  "PC16:0a/18:2", "PC16:0a/18:1", "PC18:0a/16:0", "PC16:0/18:2",
  "PC16:0/18:1", "PC18:0/16:0", "PC18:2a/18:1", "PC18:0a/18:2",

  "PC18:0a/18:1", "PC16:0/20:4", "PC18:1/18:2", "PC18:0/18:2",
  "PC18:0/18:1", "PC16:0/22:6", "PC18:1/20:4")

Phospholipids_perc_competition_LPC <- c(
  "LPC14:0", "LPC16:1", "LPC16:0", "LPC18:2", "LPC18:1",
  "LPC18:0", "LPC20:4")

Phospholipids_perc_competition_SM <- c(
  "SM16:0", "SM18:0", "SM20:0", "SM22:0", "SM24:1", "SM24:0")

Phospholipids_perc_competition_PE <- c(
  "PE16:0/18:2", "PE16:0/18:1", "PE16:1a/20:4", "PE18:1a/18:1",
  "PE18:0a/18:1", "PE16:0/20:4", "PE18:1/18:2", "PE18:1/18:1",
  "PE18:0/18:1", "PE16:0a/22:6", "PE18:1a/20:4", "PE18:0a/20:4",
  "PE16:0/22:6", "PE18:1/20:4", "PE18:0/20:4", "PE18:0/20:3",
  "PE18:0/20:2", "PE18:2a/22:6", "PE18:1/22:6", "PE18:0/22:6",
  "PE18:0/22:5", "PE18:0/22:4", "PE18:0/22:3", "PE16:1/24:1")

Phospholipids_perc_competition_PG <- c(
  "PG16:0/16:1", "PG16:0/16:0", "PG16:0a/18:1", "PG14:0/20:4",
  "PG16:1/18:2", "PG16:0/18:2", "PG16:0/18:1", "PG16:0/18:0",
  "PG18:1a/18:1", "PG18:0a/18:1", "PG18:1/18:2", "PG18:1/18:1","PG18:0/18:1")

Phospholipids_perc_competition_PI <- c(
  "PI12:0a/20:4", "PI12:0a/22:6", "PI16:0/18:2", "PI16:0/18:1",
  "PI16:0/20:4", "PI18:1/18:2", "PI18:0/18:2", "PI18:0/18:1", 
  "PI18:1/20:4", "PI18:0/20:4", "PI18:0/20:3")

Phospholipids_perc_competition_PS <- c(
  "PS16:0/18:1", "PS18:0a/18:1", "PS18:1/18:2", "PS18:0/18:2",
  "PS18:0/18:1", "PS18:0/20:4", "PS18:0/20:3", "PS20:5a/20:4", "PS20:4a/20:4", "PS18:0/22:6", "PS18:0/22:5", "PS18:0/22:4","PS20:4a/22:5", "PS20:3a/22:5"
)

# Plotting for each group of phospholipids
plot_change_from_baseline(pld_data, oxidative_stress_cols,
                          "Oxidative Stress & Inorganic Species",
                          "oxidative_stress_change_from_baseline.png")
plot_change_from_baseline(pld_data, Free_thiols_nM,
                          "Free Thiols (nM)",
                          "free_thiols_change_from_baseline.png")
plot_change_from_baseline(pld_data, Total_thiols_nM,
                          "Total Thiols (uM)",
                          "total_thiols_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_PC,
                          "Phospholipids - PC",
                          "phospholipids_PC_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_LPC,
                          "Phospholipids - LPC",
                          "phospholipids_LPC_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_SM,
                          "Phospholipids - SM",
                          "phospholipids_SM_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_PE,
                          "Phospholipids - PE",
                          "phospholipids_PE_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_PG,
                          "Phospholipids - PG",
                          "phospholipids_PG_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_PI,
                          "Phospholipids - PI",
                          "phospholipids_PI_change_from_baseline.png")
plot_change_from_baseline(pld_data, Phospholipids_perc_competition_PS,
                          "Phospholipids - PS",
                          "phospholipids_PS_change_from_baseline.png")  
