# ============================================================
# POIROT ProxiMeta ARG Analysis
# Integrates ProxiMeta ARG data with MetaPhlAn profiles
# 5 patients, paired pre/post antibiotic (T1=A, T2=B)
# Author: Eleni Kalopedis
# Date: 2026-06-10
# ============================================================

library(ggplot2)
library(dplyr)

# ============================================================
# 1. IMPORT ALL ARG FILES
# ============================================================

proxmeta_dir <- "/scratch/prj/chmi_rbiome/project/metadata/PROXIMETA_POIROT"

# Import amr_host_content — total ARG burden per MAG per sample
host_files <- list.files(proxmeta_dir, 
                          pattern="amr_host_content.tsv", 
                          full.names=TRUE)

host_list <- lapply(host_files, function(f) {
  sample <- gsub("_amr_host_content.tsv", "", basename(f))
  df <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  df$sample <- sample
  df
})
host_all <- do.call(rbind, host_list)

# Import amr_gene_to_host_summary — gene-level ARG with taxonomy
gene_files <- list.files(proxmeta_dir,
                          pattern="amr_gene_to_host_summary.tsv",
                          full.names=TRUE)

gene_list <- lapply(gene_files, function(f) {
  sample <- gsub("_amr_gene_to_host_summary.tsv", "", basename(f))
  df <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  df$sample <- sample
  df
})
gene_all <- do.call(rbind, gene_list)

# Import AMR_annotations — quality filtered ARG table
annot_files <- list.files(proxmeta_dir,
                           pattern="AMR_annotations_summary.tsv",
                           full.names=TRUE)

annot_list <- lapply(annot_files, function(f) {
  sample <- gsub("_AMR_annotations_summary.tsv", "", basename(f))
  df <- read.table(f, header=TRUE, sep="\t", 
                   stringsAsFactors=FALSE, quote="")
  df$sample <- sample
  df
})
annot_all <- do.call(rbind, annot_list)

cat("ARG host content rows:", nrow(host_all), "\n")
cat("ARG gene summary rows:", nrow(gene_all), "\n")
cat("ARG annotations rows:", nrow(annot_all), "\n")
cat("Samples:", unique(host_all$sample), "\n")

# ============================================================
# 2. DEFINE PAIRS
# ============================================================

pairs <- data.frame(
  patient = c("AKQ014", "AKQ001", "AKQ002", "AKQ003", "AKQ005"),
  T1 = c("AKQ014A", "AKQ001A", "AKQ002A", "AKQ003A", "AKQ005A"),
  T2 = c("AKQ014B", "AKQ001B", "AKQ002B", "AKQ003B", "AKQ005B"),
  stringsAsFactors=FALSE
)

# ============================================================
# 3. TOTAL ARG BURDEN PER SAMPLE
# ============================================================

arg_burden <- host_all %>%
  group_by(sample) %>%
  summarise(
    total_ARGs = sum(total_amr_genes),
    genomic_ARGs = sum(genomic_amr_genes),
    viral_ARGs = sum(viral_amr_genes),
    plasmid_ARGs = sum(plasmid_amr_genes),
    MGE_fraction = (sum(viral_amr_genes) + sum(plasmid_amr_genes)) / 
                    sum(total_amr_genes)
  )

cat("\nARG burden per sample:\n")
print(arg_burden)

# Add timepoint
arg_burden$timepoint <- ifelse(grepl("A$", arg_burden$sample),
                                "T1 Pre-antibiotic",
                                "T2 Post-antibiotic")
arg_burden$patient <- gsub("[AB]$", "", arg_burden$sample)

# ============================================================
# 4. PAIRED COMPARISON — TOTAL ARG BURDEN T1 vs T2
# ============================================================

complete_pairs <- pairs[pairs$T1 %in% arg_burden$sample & 
                         pairs$T2 %in% arg_burden$sample, ]
cat("\nComplete pairs for analysis:", nrow(complete_pairs), "\n")

if(nrow(complete_pairs) >= 3) {
  T1_burden <- arg_burden$total_ARGs[match(complete_pairs$T1, arg_burden$sample)]
  T2_burden <- arg_burden$total_ARGs[match(complete_pairs$T2, arg_burden$sample)]
  
  cat("\nT1 total ARGs:", T1_burden, "\n")
  cat("T2 total ARGs:", T2_burden, "\n")
  cat("Mean T1:", round(mean(T1_burden), 1), "\n")
  cat("Mean T2:", round(mean(T2_burden), 1), "\n")
  
  if(nrow(complete_pairs) >= 4) {
    wt <- wilcox.test(T1_burden, T2_burden, paired=TRUE, exact=FALSE)
    cat("Wilcoxon p-value:", round(wt$p.value, 4), "\n")
  }
}

# ============================================================
# 5. ARG CLASS COMPOSITION
# ============================================================

class_comp <- gene_all %>%
  group_by(sample, gene_class) %>%
  summarise(count = n(), .groups="drop")

cat("\nARG classes detected:\n")
print(sort(unique(class_comp$gene_class)))

# ============================================================
# 6. MGE-LINKED ARG BURDEN
# ============================================================

mge_burden <- gene_all %>%
  group_by(sample) %>%
  summarise(
    total = n(),
    plasmid = sum(origin == "plasmid"),
    genomic = sum(origin == "genomic"),
    viral = sum(origin == "viral"),
    mge_fraction = round((sum(origin == "plasmid") + 
                          sum(origin == "viral")) / n(), 3)
  )

cat("\nMGE-linked ARG fraction per sample:\n")
print(mge_burden)

# ============================================================
# 7. VISUALISATION — ARG BURDEN BOXPLOT
# ============================================================

p1 <- ggplot(arg_burden, 
             aes(x=timepoint, y=total_ARGs, fill=timepoint)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_point(aes(group=patient), size=3, alpha=0.8) +
  geom_line(aes(group=patient), alpha=0.4, colour="grey50") +
  scale_fill_manual(values=c("T1 Pre-antibiotic"="#3498DB",
                              "T2 Post-antibiotic"="#E74C3C")) +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Total ARG Burden - POIROT Cohort",
       subtitle="Pre vs Post-antibiotic, paired by patient",
       x="Timepoint", y="Total ARG count")

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_ARG_burden.pdf",
       plot=p1, width=6, height=6, dpi=300)

# ============================================================
# 8. VISUALISATION — ARG CLASS COMPOSITION BARPLOT
# ============================================================

class_comp$timepoint <- ifelse(grepl("A$", class_comp$sample),
                                "T1 Pre-antibiotic",
                                "T2 Post-antibiotic")
class_comp$patient <- gsub("[AB]$", "", class_comp$sample)

p2 <- ggplot(class_comp, 
             aes(x=sample, y=count, fill=gene_class)) +
  geom_bar(stat="identity") +
  facet_grid(~timepoint, scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.text=element_text(size=7),
        legend.key.size=unit(0.4, "cm"),
        plot.title=element_text(hjust=0.5, face="bold")) +
  labs(title="ARG Class Composition - POIROT Cohort",
       x="Sample", y="ARG count", fill="ARG Class")

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_ARG_classes.pdf",
       plot=p2, width=12, height=6, dpi=300)

# ============================================================
# 9. SAVE RESULTS
# ============================================================

write.csv(arg_burden,
  "/scratch/prj/chmi_rbiome/project/results/POIROT_ARG_burden.csv",
  row.names=FALSE)

write.csv(gene_all,
  "/scratch/prj/chmi_rbiome/project/results/POIROT_ARG_gene_summary.csv",
  row.names=FALSE)

cat("\nAll results saved\n")
cat("=== ARG Analysis Complete ===\n")

# ============================================================
# 10. BASELINE RESISTOME — T1 SAMPLES ONLY
# ============================================================

# Filter to T1 only
gene_T1 <- gene_all[grepl("A$", gene_all$sample), ]
annot_T1 <- annot_all[grepl("A$", annot_all$sample), ]
host_T1 <- host_all[grepl("A$", host_all$sample), ]

cat("\n=== BASELINE RESISTOME (T1 only) ===\n")
cat("Samples:", unique(gene_T1$sample), "\n")
cat("Total ARG hits:", nrow(gene_T1), "\n")
cat("Unique genes:", length(unique(gene_T1$gene_id)), "\n")

# Core resistome — ARGs present in ALL 5 T1 patients
gene_presence <- gene_T1 %>%
  group_by(gene_id) %>%
  summarise(n_patients = n_distinct(sample),
            classes = paste(unique(gene_class), collapse="/"),
            hosts = paste(unique(gsub(" \\(.*", "", bin_taxonomy)), 
                         collapse="/"),
            origin = paste(unique(origin), collapse="/"))

core <- gene_presence[gene_presence$n_patients == 5, ]
cat("\nCore resistome (present in all 5 T1 patients):", nrow(core), "ARGs\n")
print(core[, c("gene_id", "classes", "hosts", "origin")])

# Accessory — patient specific
accessory <- gene_presence[gene_presence$n_patients == 1, ]
cat("\nAccessory resistome (patient-specific):", nrow(accessory), "ARGs\n")

# ARG burden per patient at baseline
baseline_burden <- host_T1 %>%
  group_by(sample) %>%
  summarise(
    total_ARGs = sum(total_amr_genes),
    genomic_ARGs = sum(genomic_amr_genes),
    plasmid_ARGs = sum(plasmid_amr_genes),
    viral_ARGs = sum(viral_amr_genes),
    mge_fraction = round((sum(plasmid_amr_genes) + 
                          sum(viral_amr_genes)) / sum(total_amr_genes), 3)
  )

cat("\nBaseline ARG burden per patient:\n")
print(baseline_burden)

# Top ARG-carrying hosts at baseline
top_hosts <- gene_T1 %>%
  group_by(bin_taxonomy) %>%
  summarise(
    total_ARGs = n(),
    n_patients = n_distinct(sample),
    classes = paste(unique(gene_class), collapse=", ")
  ) %>%
  arrange(desc(total_ARGs))

cat("\nTop ARG-carrying hosts at baseline:\n")
print(head(top_hosts, 10))

# Visualisation — baseline ARG burden per patient
p3 <- ggplot(baseline_burden, 
             aes(x=sample, y=total_ARGs, fill=sample)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=total_ARGs), vjust=-0.3, size=3.5) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Baseline ARG Burden per Patient - POIROT Cohort",
       subtitle="Total ARGs at T1 (pre-antibiotic)",
       x="Sample", y="Total ARG count")

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_baseline_ARG_burden.pdf",
       plot=p3, width=8, height=6, dpi=300)

# Visualisation — core vs accessory resistome
prevalence_df <- data.frame(
  category = c("Core (5/5 patients)", 
               "Common (3-4 patients)",
               "Rare (2 patients)", 
               "Unique (1 patient)"),
  count = c(
    sum(gene_presence$n_patients == 5),
    sum(gene_presence$n_patients %in% 3:4),
    sum(gene_presence$n_patients == 2),
    sum(gene_presence$n_patients == 1)
  )
)

p4 <- ggplot(prevalence_df, aes(x=category, y=count, fill=category)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=count), vjust=-0.3, size=4) +
  scale_fill_manual(values=c(
    "Core (5/5 patients)"="#2ECC71",
    "Common (3-4 patients)"="#3498DB",
    "Rare (2 patients)"="#E67E22",
    "Unique (1 patient)"="#E74C3C"
  )) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=30, hjust=1),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Core vs Accessory Resistome - POIROT Cohort",
       subtitle="ARG prevalence across 5 patients at baseline (T1)",
       x="", y="Number of ARGs")

ggsave("/scratch/prj/chmi_rbiome/project/figures/POIROT_core_accessory_resistome.pdf",
       plot=p4, width=8, height=6, dpi=300)

# Save baseline results
write.csv(gene_presence,
  "/scratch/prj/chmi_rbiome/project/results/POIROT_baseline_resistome.csv",
  row.names=FALSE)

write.csv(core,
  "/scratch/prj/chmi_rbiome/project/results/POIROT_core_resistome.csv",
  row.names=FALSE)

cat("\nBaseline resistome figures and results saved\n")
