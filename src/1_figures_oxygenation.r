# produce figures for oxygenation data

library(tidyverse)
library(patchwork)

# run cleaning script if processed data does not exist
if (!file.exists("data/processed/cleaned_oxygenation_data.rds")) {
    message("Cleaned data not found — running src/0_clean_oxygenation.r ...")
    dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
    source("src/0_clean_oxygenation.r")
}

o2_data_long <- readRDS("data/processed/cleaned_oxygenation_data.rds")

# --- helper: produce a Lowess/LM split plot with risk table -------------------
plot_oxygenation <- function(data, var_name, y_label, y_breaks, filename) {

    dat <- data %>% filter(variable == var_name)
    n_treatments <- n_distinct(dat$Treatment)
    treatment_levels <- levels(factor(dat$Treatment))
    treatment_cols <- setNames(
        RColorBrewer::brewer.pal(max(3, n_treatments), "Set1")[seq_len(n_treatments)],
        treatment_levels
    )

    x_breaks <- seq(-12, 72, by = 12)
    x_limits <- c(-14, 78)

    # --- main plot ---
    p <- ggplot(dat, aes(x = TimeSinceRandomisation, y = value)) +
        geom_smooth(
            data = dat %>% filter(TimeSinceRandomisation < 0),
            aes(colour = "Pre-randomisation"),
            method = "lm", se = TRUE, fill = "grey70"
        ) +
        geom_smooth(
            data = dat %>% filter(TimeSinceRandomisation >= 0),
            aes(colour = factor(Treatment)),
            method = "loess", se = TRUE, span = 0.75, fill = "grey70"
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
        scale_colour_manual(
            name = "Treatment group",
            values = c("Pre-randomisation" = "grey40", treatment_cols)
        ) +
        scale_x_continuous(
            limits = x_limits,
            breaks = x_breaks,
            labels = x_breaks,
            expand = expansion(mult = c(0, 0))
        ) +
        scale_y_continuous(
            breaks = y_breaks,
            expand = expansion(mult = c(0, 0))
        ) +
        coord_cartesian(ylim = range(y_breaks)) +
        labs(x = NULL, y = y_label) +
        theme_bw(base_size = 16) +
        theme(
            legend.position = "inside", legend.position.inside = c(0.97, 0.97),
            legend.justification = c(1, 1),
            legend.background = element_rect(fill = alpha("white", 0)),
            legend.key = element_rect(fill = NA),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 16),
            axis.title.y = element_text(margin = margin(r = 0)),
            axis.text = element_text(size = 14),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(5.5, 5.5, 0, 5.5)
        )

    # --- risk table (exclude time = -12) ---
    time_breaks <- seq(0, 72, by = 12)
    treatment_levels_sorted <- sort(unique(dat$Treatment))
    risk <- expand_grid(Treatment = treatment_levels_sorted, time = time_breaks) %>%
        mutate(n = map2_dbl(Treatment, time, function(trt, t) {
            n_distinct(dat$MECROXStudy[dat$Treatment == trt & dat$TimeSinceRandomisation >= t])
        }))

    p_risk <- ggplot(risk, aes(x = time, y = Treatment, label = as.character(n))) +
        geom_text(size = 4.5, fontface = "bold") +
        scale_x_continuous(
            limits = x_limits,
            breaks = x_breaks,
            labels = x_breaks,
            expand = expansion(mult = c(0, 0))
        ) +
        labs(x = "Time since randomisation (hours)", y = NULL) +
        theme_bw(base_size = 16) +
        theme(
            axis.title.x = element_text(size = 16),
            axis.title.y = element_blank(),
            axis.text = element_text(size = 14),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(0, 5.5, 5.5, 5.5)
        )

    # --- combine (no axes collect — align via matching margins) ---
    combined <- p / p_risk + plot_layout(heights = c(5, 1))

    ggsave(paste0("output/figures/", filename), plot = combined, width = 8, height = 6.5, dpi = 300)
    message("Saved: ", filename)
}

# --- generate figures ----------------------------------------------------------

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

plot_oxygenation(o2_data_long, "FiO2",    "FiO2 (%)",     seq(30, 70, by = 10), "lowess_fio2_split.png")
plot_oxygenation(o2_data_long, "PO2",     "PaO2 (kPa)",   seq(8, 16, by = 2),   "lowess_pao2_split.png")
plot_oxygenation(o2_data_long, "PFRatio", "P/F ratio",    seq(20, 32, by = 2),  "lowess_pfratio_split.png")
plot_oxygenation(o2_data_long, "SpO2",    "SpO2 (%)",     seq(88, 100, by = 2), "lowess_spo2_split.png")
