# ============================================================
# PROFIT Baseline Resistome Analysis
# ProxiMeta Explorer amr_composition.tsv files
# 17 patients, Day 1 baseline
# Author: Eleni Kalopedis
# Date: 2026-06-13
# ============================================================

library(dplyr)
library(ggplot2)

# ============================================================
# 1. IMPORT ALL 17 PATIENT AMR FILES
# ============================================================
base <- "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_ProxiMeta/PME_PROFIT"
folders <- list.dirs(base, recursive=FALSE)

amr_list <- lapply(folders, function(f) {
  patient <- basename(f)
  amr_file <- file.path(f, "amr_composition.tsv")
  if(file.exists(amr_file)) {
    df <- read.table(amr_file, header=TRUE, sep="\t",
                     stringsAsFactors=FALSE, quote="")
    df$patient_folder <- patient
    df$patient_id <- gsub("(p\\d+)-d1.*", "\\1", patient)
    df
  }
})

amr_all <- do.call(rbind, amr_list[!sapply(amr_list, is.null)])
amr_only <- amr_all[amr_all$Element.type == "AMR",]

cat("Total rows:", nrow(amr_all), "\n")
cat("Patients:", length(unique(amr_all$patient_id)), "\n")
cat("Element types:\n")
print(table(amr_all$Element.type))

# ============================================================
# 2. ARG BURDEN PER PATIENT
# ============================================================
burden <- amr_only %>%
  group_by(patient_id) %>%
  summarise(
    total_ARG = n(),
    plasmid_ARGs = sum(Gene.contig.type == "plasmid", na.rm=TRUE),
    genomic_ARGs = sum(Gene.contig.type == "genomic", na.rm=TRUE),
    mge_fraction = round(plasmid_ARGs / total_ARG * 100, 1),
    unique_genes = n_distinct(Gene.ID),
    .groups="drop"
  ) %>%
  arrange(desc(total_ARG))

cat("\nPROFIT baseline ARG burden:\n")
cat("Mean:", round(mean(burden$total_ARG),1), "\n")
cat("Range:", min(burden$total_ARG), "-", max(burden$total_ARG), "\n")
cat("Mean MGE fraction:", round(mean(burden$mge_fraction),1), "%\n")

# ============================================================
# 3. CROSS-COHORT COMPARISON — PROFIT vs POIROT
# ============================================================
resist_poirot <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROXIMETA_POIROT/POIROTallsamplesresistance.csv",
  stringsAsFactors=FALSE, check.names=FALSE)

T1_samples <- c("AKQ_005_A","AKQ_003_A","AKQ_014_A","AKQ_002_A")
poirot_pmg <- resist_poirot[resist_poirot$`Sample ID` %in% T1_samples &
                              resist_poirot$`Element type` == "AMR",]

poirot_burden <- poirot_pmg %>%
  group_by(`Sample ID`) %>%
  summarise(total_ARG=n(),
            plasmid=sum(`Gene contig type`=="plasmid"),
            mge_fraction=round(plasmid/n()*100,1),
            .groups="drop")

wt <- wilcox.test(burden$total_ARG, poirot_burden$total_ARG, exact=FALSE)
cat("\nWilcoxon PROFIT vs POIROT ARG burden p-value:", round(wt$p.value,4), "\n")

# ============================================================
# 4. CORE VS ACCESSORY RESISTOME
# ============================================================
gene_presence <- amr_only %>%
  group_by(Gene.ID) %>%
  summarise(n_patients=n_distinct(patient_id),
            classes=paste(unique(Gene.Class), collapse="/"),
            .groups="drop")

prevalence_df <- data.frame(
  category=c("Core (17/17)","Common (10-16)","Rare (2-9)","Unique (1/17)"),
  count=c(
    sum(gene_presence$n_patients == 17),
    sum(gene_presence$n_patients >= 10 & gene_presence$n_patients < 17),
    sum(gene_presence$n_patients >= 2 & gene_presence$n_patients < 10),
    sum(gene_presence$n_patients == 1)
  )
)

core_gene <- gene_presence[gene_presence$n_patients == 17,]
cat("\nCore gene (all 17 patients):", core_gene$Gene.ID, "\n")

# ============================================================
# 5. TOP ARG-CARRYING SPECIES
# ============================================================
top_species <- amr_only %>%
  filter(!grepl("NoHost|Unclassified", Species) & Species != "") %>%
  group_by(Species) %>%
  summarise(total_ARGs=n(), n_patients=n_distinct(patient_id),
            .groups="drop") %>%
  arrange(desc(total_ARGs))

# ============================================================
# 6. FIGURES
# ============================================================
plot_df <- rbind(
  data.frame(cohort="PROFIT\n(chronic, n=17)", total_ARG=burden$total_ARG),
  data.frame(cohort="POIROT\n(acute, n=4)", total_ARG=poirot_burden$total_ARG)
)

p1 <- ggplot(plot_df, aes(x=cohort, y=total_ARG, fill=cohort)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(width=0.15, size=2.5, alpha=0.7) +
  scale_fill_manual(values=c("PROFIT\n(chronic, n=17)"="#E74C3C",
                              "POIROT\n(acute, n=4)"="#3498DB")) +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Baseline ARG Burden — Cross-Cohort Comparison",
       subtitle="Wilcoxon p=1.0 (ns)",
       x="Cohort", y="Total AMR gene count")

ggsave("/scratch/prj/chmi_rbiome/project/figures/CrossCohort_ARG_burden.pdf",
       p1, width=6, height=6, dpi=300)

cat("\nAnalysis complete\n")
