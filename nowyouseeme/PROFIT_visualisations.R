# ============================================================
# PROFIT Clinical Analysis — Visualisations
# Reads from saved correlation results CSV
# Author: Eleni Kalopedis
# Date: 2026-05-31
# ============================================================

library(ggplot2)
library(dplyr)

# ============================================================
# LOAD DATA
# ============================================================

# Load correlation results
results <- read.csv(
  "/scratch/prj/chmi_rbiome/project/results/PROFIT_correlations_all.csv",
  stringsAsFactors = FALSE
)

# Load baseline data for scatter plots
baseline <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Baseline.csv.csv",
  skip = 2, header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)

# Rename key columns
colnames(baseline)[1]  <- "patient_id"
colnames(baseline)[2]  <- "pseudonymised_id"
colnames(baseline)[9]  <- "MELD"
colnames(baseline)[52] <- "faecal_IL22"
colnames(baseline)[58] <- "faecal_IL1b"
colnames(baseline)[60] <- "faecal_IL6"
colnames(baseline)[78] <- "faecal_ammonia"
colnames(baseline)[42] <- "Calprotectin"
colnames(baseline)[82] <- "E_faecalis_copies"
colnames(baseline)[48] <- "faecal_IL17F"

# Force numeric
for(col in c("faecal_IL22", "faecal_IL1b", "faecal_IL6",
             "faecal_ammonia", "Calprotectin",
             "E_faecalis_copies", "faecal_IL17F", "MELD")) {
  baseline[[col]] <- as.numeric(as.character(baseline[[col]]))
}

# Filter arms
fmt     <- baseline[which(baseline$IMP == 1), ]
placebo <- baseline[which(baseline$IMP == 2), ]

# ============================================================
# FIGURE 1 — CORRELATION HEATMAP
# Nominally significant novel correlations
# ============================================================

exploratory <- results[!is.na(results$p_value) &
                         results$p_value <= 0.05, ]
exploratory$arm_tp <- paste(exploratory$arm,
                             exploratory$timepoint)

p1 <- ggplot(exploratory,
             aes(x = cytokine, y = outcome, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(rho, "\np=", p_value)),
            size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white",
                       high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1),
                       name = "Spearman\nrho") +
  facet_wrap(~arm_tp, ncol = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(title = "Novel Clinical Correlations — PROFIT Trial",
       subtitle = "Nominally significant associations (p < 0.05) not reported in pre-print",
       x = "Faecal Cytokine", y = "Outcome Variable")

ggsave("/scratch/prj/chmi_rbiome/project/figures/Fig1_PROFIT_correlation_heatmap.pdf",
       plot = p1, width = 14, height = 10, dpi = 300)

cat("Figure 1 saved\n")

# ============================================================
# FIGURE 2 — FAECAL AMMONIA vs IL-1b ACROSS TIMEPOINTS
# FMT arm — novel finding at D90
# ============================================================

ammonia_il1b <- fmt[, c("pseudonymised_id", "faecal_ammonia",
                         "faecal_IL1b")]
ammonia_il1b$timepoint <- factor(
  gsub(".*_", "", fmt$pseudonymised_id),
  levels = c("D0", "D7", "D30", "D90")
)
ammonia_il1b <- na.omit(ammonia_il1b)

p2 <- ggplot(ammonia_il1b,
             aes(x = faecal_ammonia, y = faecal_IL1b)) +
  geom_point(colour = "#E74C3C", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#2C3E50", linewidth = 0.8,
              fill = "#BDC3C7") +
  facet_wrap(~timepoint, ncol = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        strip.text = element_text(face = "bold")) +
  labs(title = "Faecal Ammonia vs Faecal IL-1β — FMT Group",
       subtitle = "Spearman rho = -0.736 at D90 (p = 0.002)",
       x = "Faecal Ammonia (nMol)",
       y = "Faecal IL-1β (normalised to protein)")

ggsave("/scratch/prj/chmi_rbiome/project/figures/Fig2_PROFIT_ammonia_IL1b_FMT.pdf",
       plot = p2, width = 12, height = 4, dpi = 300)

cat("Figure 2 saved\n")

# ============================================================
# FIGURE 3 — CALPROTECTIN vs IL-6 AT D30
# FMT arm — novel finding
# ============================================================

calp_il6 <- fmt[grep("_D30$", fmt$pseudonymised_id),
                c("Calprotectin", "faecal_IL6")]
calp_il6 <- na.omit(calp_il6)

p3 <- ggplot(calp_il6,
             aes(x = Calprotectin, y = faecal_IL6)) +
  geom_point(colour = "#3498DB", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#2C3E50", linewidth = 0.8,
              fill = "#BDC3C7") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(title = "Calprotectin vs Faecal IL-6 at D30 — FMT Group",
       subtitle = "Spearman rho = 0.529, p = 0.043, n = 15",
       x = "Faecal Calprotectin (µg/g)",
       y = "Faecal IL-6 (normalised to protein)")

ggsave("/scratch/prj/chmi_rbiome/project/figures/Fig3_PROFIT_calprotectin_IL6_D30.pdf",
       plot = p3, width = 7, height = 6, dpi = 300)

cat("Figure 3 saved\n")

# ============================================================
# FIGURE 4 — E. FAECALIS vs IL-17F AT D0
# Placebo arm — FDR significant finding
# ============================================================

efaec_il17f <- placebo[grep("_D0$", placebo$pseudonymised_id),
                        c("E_faecalis_copies", "faecal_IL17F")]
efaec_il17f <- na.omit(efaec_il17f)

p4 <- ggplot(efaec_il17f,
             aes(x = E_faecalis_copies, y = faecal_IL17F)) +
  geom_point(colour = "#9B59B6", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#2C3E50", linewidth = 0.8,
              fill = "#BDC3C7") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(title = "E. faecalis Copy Number vs Faecal IL-17F at Baseline",
       subtitle = "Placebo group — Spearman rho = 1.0, p_adj < 0.001, n = 5",
       x = "E. faecalis DNA Copy Number",
       y = "Faecal IL-17F (normalised to protein)")

ggsave("/scratch/prj/chmi_rbiome/project/figures/Fig4_PROFIT_Efaecalis_IL17F_D0.pdf",
       plot = p4, width = 7, height = 6, dpi = 300)

cat("Figure 4 saved\n")
cat("\n=== All figures saved to figures folder ===\n")

