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

# ============================================================
# 8. DIFFERENTIAL ABUNDANCE — WILCOXON SIGNED-RANK
# ============================================================

results_da <- data.frame(
  species = character(), mean_T1 = numeric(),
  mean_T2 = numeric(), log2FC = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
)

for(i in 1:nrow(species)) {
  T1_vals <- as.numeric(species[i, pairs$T1])
  T2_vals <- as.numeric(species[i, pairs$T2])
  
  if(max(c(T1_vals, T2_vals)) > 0.1) {
    test <- wilcox.test(T1_vals, T2_vals, paired=TRUE, exact=FALSE)
    mean_T1 <- mean(T1_vals)
    mean_T2 <- mean(T2_vals)
    lfc <- log2((mean_T2 + 0.001) / (mean_T1 + 0.001))
    
    results_da <- rbind(results_da, data.frame(
      species = species$clade_name[i],
      mean_T1 = round(mean_T1, 4),
      mean_T2 = round(mean_T2, 4),
      log2FC = round(lfc, 3),
      p_value = round(test$p.value, 4),
      stringsAsFactors = FALSE
    ))
  }
}

results_da$p_adjusted <- round(p.adjust(results_da$p_value, method="BH"), 4)
results_da$species_short <- gsub(".*s__", "", results_da$species)

# Note: minimum achievable p-value with n=5 pairs is 0.0625
# Statistical significance is mathematically unachievable at p<0.05
cat("Minimum p-value achieved:", min(results_da$p_value), "\n")
cat("Species tested:", nrow(results_da), "\n")

# Save results
write.csv(results_da,
  "/scratch/prj/chmi_rbiome/project/results/POIROT_differential_abundance.csv",
  row.names=FALSE)

# Volcano plot
results_da$colour <- ifelse(results_da$log2FC > 1 & results_da$p_value < 0.1,
                             "Increased",
                      ifelse(results_da$log2FC < -1 & results_da$p_value < 0.1,
                             "Decreased", "Not significant"))

p3 <- ggplot(results_da, aes(x=log2FC, y=-log10(p_value),
                              colour=colour)) +
  geom_point(alpha=0.6, size=2) +
  geom_text(data=head(results_da[order(results_da$p_value),], 8),
            aes(label=species_short), size=2.5, hjust=-0.1) +
  scale_colour_manual(values=c("Increased"="#E74C3C",
                                "Decreased"="#3498DB",
                                "Not significant"="grey60")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", alpha=0.5) +
  theme_bw() +
  theme(legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Differential Abundance - POIROT Cohort",
       subtitle="Post vs Pre-antibiotic (T2 vs T1), n=5 pairs",
       x="log2 Fold Change (T2/T1)",
       y="-log10(p-value)")

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_differential_abundance.pdf",
       plot=p3, width=10, height=7, dpi=300)

cat("Volcano plot saved\n")

# ============================================================
# 9. PERMANOVA — TIMEPOINT AND PATIENT
# ============================================================

permanova_manual <- function(dist_mat, groups, n_perm=999) {
  calc_F <- function(d, g) {
    groups_unique <- unique(g)
    n <- length(g)
    SS_total <- sum(d^2) / n
    SS_within <- 0
    for(grp in groups_unique) {
      idx <- which(g == grp)
      n_grp <- length(idx)
      if(n_grp > 1) SS_within <- SS_within + sum(d[idx,idx]^2) / n_grp
    }
    SS_between <- SS_total - SS_within
    df_between <- length(groups_unique) - 1
    df_within <- n - length(groups_unique)
    (SS_between/df_between) / (SS_within/df_within)
  }
  F_obs <- calc_F(dist_mat, groups)
  F_perm <- replicate(n_perm, calc_F(dist_mat, sample(groups)))
  p_val <- (sum(F_perm >= F_obs) + 1) / (n_perm + 1)
  list(F=round(F_obs,3), p=round(p_val,4))
}

timepoint_groups <- ifelse(grepl("A$", rownames(dist_mat)), "T1", "T2")
result_tp <- permanova_manual(dist_mat, timepoint_groups)
cat("PERMANOVA by timepoint - F:", result_tp$F, "p:", result_tp$p, "\n")

patient_groups <- gsub("KINGCO-007-", "", gsub("[AB]$", "", rownames(dist_mat)))
result_pt <- permanova_manual(dist_mat, patient_groups)
cat("PERMANOVA by patient - F:", result_pt$F, "p:", result_pt$p, "\n")

permanova_results <- data.frame(
  test = c("Timepoint", "Patient"),
  F_statistic = c(result_tp$F, result_pt$F),
  p_value = c(result_tp$p, result_pt$p)
)

write.csv(permanova_results,
  "/scratch/prj/chmi_rbiome/project/results/POIROT_PERMANOVA_results.csv",
  row.names=FALSE)

cat("PERMANOVA results saved\n")

# ============================================================
# 10. TOP SPECIES ABUNDANCE BARPLOT
# ============================================================

species$mean_abund <- rowMeans(species[,2:11])
top15 <- species[order(species$mean_abund, decreasing=TRUE), ][1:15, ]
top15$species_short <- gsub(".*s__", "", top15$clade_name)

long_df <- data.frame()
for(i in 1:nrow(top15)) {
  for(j in 2:11) {
    long_df <- rbind(long_df, data.frame(
      species = top15$species_short[i],
      sample = colnames(top15)[j],
      abundance = as.numeric(top15[i,j]),
      stringsAsFactors = FALSE
    ))
  }
}

long_df$timepoint <- ifelse(grepl("A$", long_df$sample),
                             "T1 Pre-antibiotic", "T2 Post-antibiotic")
long_df$patient <- gsub("KINGCO-007-", "", gsub("[AB]$", "", long_df$sample))
long_df$sample_label <- paste0(long_df$patient, "\n", long_df$timepoint)

species_order <- top15$species_short[order(top15$mean_abund, decreasing=TRUE)]
long_df$species <- factor(long_df$species, levels=rev(species_order))

p4 <- ggplot(long_df, aes(x=sample_label, y=abundance, fill=species)) +
  geom_bar(stat="identity") +
  facet_grid(~timepoint, scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        legend.text=element_text(size=7),
        legend.key.size=unit(0.4, "cm"),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Top 15 Species Abundance - POIROT Cohort",
       subtitle="Pre vs Post-antibiotic",
       x="Sample", y="Relative Abundance (%)",
       fill="Species")

# ============================================================
# OBSERVED RICHNESS + SHANNON COMBINED FIGURE
# Postdoc requested observed diversity alongside Shannon
# Date: 2026-06-20
# ============================================================
# Key finding: AKQ005 is the unique outlier for OBSERVED richness
# (423->281, -142 species) consistent with established mass
# extinction phenotype (174 extinctions, E.coli/T4SS elimination)
# AKQ001 is the unique outlier for SHANNON (4.12->3.82 decrease)
# Divergence between metrics demonstrates richness and evenness
# capture different ecological dimensions
# Observed paired Wilcoxon: p=0.4185 (NS)
# Shannon paired Wilcoxon: p=0.2807 (NS)

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_species_barplot.pdf",
       plot=p4, width=14, height=8, dpi=300)

cat("Species barplot saved\n")

# ============================================================
# 11. ADDITIONAL ANALYSES — CRITICAL EVALUATION
# ============================================================

# Shannon components — richness and evenness
richness <- function(x) sum(x > 0)
evenness <- function(x) { h<-shannon(x); s<-richness(x); if(s>1) h/log(s) else 0 }

richness_scores <- sapply(species_mat, richness)
evenness_scores <- sapply(species_mat, evenness)

div_df <- data.frame(
  sample=names(shannon_scores),
  shannon=round(shannon_scores,3),
  richness=richness_scores,
  evenness=round(evenness_scores,3),
  timepoint=ifelse(grepl("A$",names(shannon_scores)),"T1","T2"),
  patient=gsub("KINGCO-007-","",gsub("[AB]$","",names(shannon_scores)))
)

cat("T1 means - Shannon:", round(mean(div_df$shannon[div_df$timepoint=="T1"]),3),
    "| Richness:", round(mean(div_df$richness[div_df$timepoint=="T1"]),1),
    "| Evenness:", round(mean(div_df$evenness[div_df$timepoint=="T1"]),3),"\n")
cat("T2 means - Shannon:", round(mean(div_df$shannon[div_df$timepoint=="T2"]),3),
    "| Richness:", round(mean(div_df$richness[div_df$timepoint=="T2"]),1),
    "| Evenness:", round(mean(div_df$evenness[div_df$timepoint=="T2"]),3),"\n")

# Unclassified fraction vs Shannon correlation
unclassified <- poirot[poirot$clade_name=="UNCLASSIFIED",]
unclassified_vec <- as.numeric(unclassified[,2:11])
names(unclassified_vec) <- colnames(unclassified)[2:11]
cor_result <- cor.test(unclassified_vec, shannon_scores, method="spearman", exact=FALSE)
cat("\nUnclassified fraction vs Shannon - rho:",
    round(cor_result$estimate,3), "p:", round(cor_result$p.value,4),"\n")

# Within vs between patient distance ratio
within_dists <- sapply(1:nrow(pairs), function(i) dist_mat[pairs$T1[i],pairs$T2[i]])
between_dists <- c()
for(i in 1:nrow(pairs)) for(j in 1:nrow(pairs)) {
  if(i!=j) between_dists <- c(between_dists, dist_mat[pairs$T1[i],pairs$T1[j]])
}
cat("\nMean within-patient distance:", round(mean(within_dists),3),"\n")
cat("Mean between-patient distance:", round(mean(between_dists),3),"\n")
cat("Ratio (between/within):", round(mean(between_dists)/mean(within_dists),2),"x\n")

# ESKAPE pathogen check
pathogens <- c("Clostridioides_difficile","Klebsiella_pneumoniae",
               "Enterococcus_faecium","Staphylococcus_aureus",
               "Acinetobacter_baumannii","Pseudomonas_aeruginosa")
cat("\n=== ESKAPE pathogen check ===\n")
for(p in pathogens) {
  hits <- species[grep(p, species$clade_name, ignore.case=TRUE),]
  if(nrow(hits)>0) {
    cat(p,"FOUND - T1 samples:",
        sum(as.numeric(hits[,pairs$T1])>0),
        "| T2 samples:",
        sum(as.numeric(hits[,pairs$T2])>0),"\n")
  } else cat(p,": not detected\n")
}

# Colonisers and extinctions
cat("\n=== Colonisers and extinctions ===\n")
for(i in 1:nrow(pairs)) {
  T1_present <- species$clade_name[species[,pairs$T1[i]]>0]
  T2_present <- species$clade_name[species[,pairs$T2[i]]>0]
  cat(pairs$patient[i],
      "- Extinctions:", length(setdiff(T1_present,T2_present)),
      "| New colonisers:", length(setdiff(T2_present,T1_present)),"\n")
}

# ============================================================
# FINAL ESKAPE FISHER TEST — Complete n=5 paired dataset
# Includes AKQ001B and AKQ002B (added 2026-06-20)
# CORRECTED: previous note estimated p≈0.047 — this was
# speculative and incorrect. Actual result below.
# ============================================================
# T1: 0/5 ESKAPE+ | T2: 3/5 ESKAPE+ (AKQ03B, AKQ05B, AKQ014B)
# Fisher exact test p-value: 0.1667 (NS)
# Report as directional trend, not significant finding,
# due to n=5 power constraint

