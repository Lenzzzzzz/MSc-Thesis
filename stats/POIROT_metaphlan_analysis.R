# ============================================================
# POIROT MetaPhlAn Analysis
# Taxonomic profiling of pre/post antibiotic samples
# KINGCO-007 cohort — 5 patients, 10 samples (A=T1, B=T2)
# Author: Eleni Kalopedis
# Date: 2026-06-10
# ============================================================

library(ggplot2)
library(dplyr)

# ============================================================
# 1. IMPORT MERGED METAPHLAN TABLE
# ============================================================

poirot <- read.table(
  "/scratch/prj/chmi_rbiome/project/results/POIROT_merged_profiles.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Dimensions:", dim(poirot), "\n")

# ============================================================
# 2. FILTER TO SPECIES LEVEL
# ============================================================

species <- poirot[grep("s__", poirot$clade_name), ]
species <- species[!grepl("t__", species$clade_name), ]
cat("Number of species:", nrow(species), "\n")

# ============================================================
# 3. CHECK UNCLASSIFIED FRACTION
# ============================================================

unclassified <- poirot[poirot$clade_name == "UNCLASSIFIED", ]
cat("\nUnclassified fraction per sample:\n")
print(round(unclassified[, 2:11], 2))

# ============================================================
# 4. DEFINE PATIENT PAIRS
# ============================================================

pairs <- data.frame(
  patient = c("AKQ014", "AKQ01", "AKQ02", "AKQ03", "AKQ05"),
  T1 = c("KINGCO-007-AKQ014A", "KINGCO-007-AKQ01A",
         "KINGCO-007-AKQ02A", "KINGCO-007-AKQ03A", "KINGCO-007-AKQ05A"),
  T2 = c("KINGCO-007-AKQ014B", "KINGCO-007-AKQ01B",
         "KINGCO-007-AKQ02B", "KINGCO-007-AKQ03B", "KINGCO-007-AKQ05B"),
  stringsAsFactors = FALSE
)

# ============================================================
# 5. ALPHA DIVERSITY — SHANNON INDEX
# ============================================================

shannon <- function(x) {
  x <- x[x > 0]
  x <- x / sum(x)
  -sum(x * log(x))
}

species_mat <- species[, 2:11]
rownames(species_mat) <- species$clade_name

shannon_scores <- sapply(species_mat, shannon)
cat("\nShannon diversity per sample:\n")
print(round(shannon_scores, 3))

T1_shannon <- shannon_scores[pairs$T1]
T2_shannon <- shannon_scores[pairs$T2]

cat("\nT1 mean Shannon:", round(mean(T1_shannon), 3), "\n")
cat("T2 mean Shannon:", round(mean(T2_shannon), 3), "\n")

wilcox_result <- wilcox.test(T1_shannon, T2_shannon, paired = TRUE)
cat("Wilcoxon p-value:", round(wilcox_result$p.value, 4), "\n")

# ============================================================
# 6. VISUALISATION — SHANNON DIVERSITY BOXPLOT
# ============================================================

shannon_df <- data.frame(
  sample = names(shannon_scores),
  shannon = shannon_scores,
  timepoint = ifelse(grepl("A$", names(shannon_scores)), "T1 Pre-antibiotic", "T2 Post-antibiotic"),
  patient = gsub("KINGCO-007-AKQ0?", "AKQ", 
                 gsub("[AB]$", "", names(shannon_scores)))
)

p1 <- ggplot(shannon_df, aes(x = timepoint, y = shannon, fill = timepoint)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(group = patient), size = 3, alpha = 0.8) +
  geom_line(aes(group = patient), alpha = 0.4, colour = "grey50") +
  scale_fill_manual(values = c("T1 Pre-antibiotic" = "#3498DB",
                               "T2 Post-antibiotic" = "#E74C3C")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  labs(title = "Alpha Diversity — POIROT Cohort",
       subtitle = paste0("Wilcoxon signed-rank p = ",
                         round(wilcox_result$p.value, 4), ", n = 5 pairs"),
       x = "Timepoint",
       y = "Shannon Diversity Index")

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_shannon_diversity.pdf",
       plot = p1, width = 6, height = 6, dpi = 300)

cat("\nShannon diversity figure saved\n")
cat("\n=== Analysis complete ===\n")

# ============================================================
# 7. BETA DIVERSITY — BRAY-CURTIS + PCoA
# ============================================================

species_mat <- t(as.matrix(species[, 2:11]))
species_mat <- species_mat / rowSums(species_mat)

bray_curtis <- function(x, y) {
  sum(abs(x - y)) / sum(x + y)
}

n <- nrow(species_mat)
dist_mat <- matrix(0, n, n)
rownames(dist_mat) <- colnames(dist_mat) <- rownames(species_mat)

for(i in 1:n) {
  for(j in 1:n) {
    dist_mat[i,j] <- bray_curtis(species_mat[i,], species_mat[j,])
  }
}

pcoa <- cmdscale(dist_mat, k=2, eig=TRUE)
var_explained <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

pcoa_df <- data.frame(
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2],
  sample = rownames(pcoa$points),
  timepoint = ifelse(grepl("A$", rownames(pcoa$points)),
                     "T1 Pre-antibiotic", "T2 Post-antibiotic"),
  patient = gsub("KINGCO-007-", "", gsub("[AB]$", "", rownames(pcoa$points)))
)

cat("Variance explained PC1:", var_explained[1], "%\n")
cat("Variance explained PC2:", var_explained[2], "%\n")

p2 <- ggplot(pcoa_df, aes(x=PC1, y=PC2, colour=timepoint, label=patient)) +
  geom_point(size=4, alpha=0.8) +
  geom_line(aes(group=patient), colour="grey50", alpha=0.5) +
  geom_text(nudge_y=0.02, size=3) +
  scale_colour_manual(values=c("T1 Pre-antibiotic"="#3498DB",
                               "T2 Post-antibiotic"="#E74C3C")) +
  theme_bw() +
  theme(legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Beta Diversity PCoA - POIROT Cohort",
       subtitle="Bray-Curtis dissimilarity, lines connect paired samples",
       x=paste0("PC1 (", var_explained[1], "% variance)"),
       y=paste0("PC2 (", var_explained[2], "% variance)"))

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_beta_diversity_PCoA.pdf",
       plot=p2, width=8, height=6, dpi=300)

cat("PCoA figure saved\n")
