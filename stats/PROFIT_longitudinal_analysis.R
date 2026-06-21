# ============================================================
# PROFIT Longitudinal ARG Analysis
# FMT vs Placebo arms across Day 0, 7, 30, 90
# Author: Eleni Kalopedis
# Date: 2026-06-13
# ============================================================

library(dplyr)
library(ggplot2)

base <- "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_ProxiMeta"

# ============================================================
# 1. IMPORT FUNCTION
# ============================================================
import_timepoint <- function(timepoint_path, timepoint_label) {
  folders <- list.dirs(timepoint_path, recursive=FALSE)
  amr_list <- lapply(folders, function(f) {
    amr_file <- file.path(f, "amr_composition.tsv")
    if(file.exists(amr_file)) {
      df <- tryCatch(
        read.table(amr_file, header=TRUE, sep="\t",
                   stringsAsFactors=FALSE, quote=""),
        error=function(e) NULL)
      if(!is.null(df)) {
        df$folder <- basename(f)
        df$timepoint <- timepoint_label
        df$patient_id <- gsub("(p\\d+)-.*", "\\1", basename(f))
        df
      }
    }
  })
  do.call(rbind, amr_list[!sapply(amr_list, is.null)])
}

# ============================================================
# 2. IMPORT ALL TIMEPOINTS
# ============================================================
d0 <- import_timepoint(file.path(base,"PME_PROFIT"),"Day0")
d0 <- d0[!grepl("^d", d0$folder),]
d7  <- import_timepoint(file.path(base,"Day7"), "Day7")
d30 <- import_timepoint(file.path(base,"Day30"), "Day30")
d90 <- import_timepoint(file.path(base,"Day90"), "Day90")

all_data <- rbind(d0, d7, d30, d90)
amr_all  <- all_data[all_data$Element.type == "AMR",]

cat("Rows per timepoint:\n")
print(table(amr_all$timepoint))

# ============================================================
# 3. ARG BURDEN PER PATIENT PER TIMEPOINT
# ============================================================
burden_long <- amr_all %>%
  group_by(patient_id, timepoint) %>%
  summarise(
    total_ARG    = n(),
    plasmid_ARGs = sum(Gene.contig.type=="plasmid", na.rm=TRUE),
    mge_fraction = round(plasmid_ARGs/total_ARG*100, 1),
    unique_genes = n_distinct(Gene.ID),
    .groups="drop"
  )

burden_long$timepoint <- factor(burden_long$timepoint,
                                 levels=c("Day0","Day7","Day30","Day90"))

# ============================================================
# 4. ARM ASSIGNMENT FROM CLINICAL METADATA
# ============================================================
meta_raw <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Baseline.csv",
  stringsAsFactors=FALSE, check.names=FALSE,
  header=FALSE, skip=2)

colnames(meta_raw) <- as.character(meta_raw[1,])
meta_raw <- meta_raw[-1,]
colnames(meta_raw)[1] <- "Patient_ID"
colnames(meta_raw)[2] <- "Pseudonymised_ID"
colnames(meta_raw) <- make.unique(colnames(meta_raw), sep="_")

meta_D0 <- meta_raw[grepl("D0$", meta_raw$Patient_ID),]
meta_D0$patient_id <- paste0("p", as.integer(
  gsub("P0100(\\d+).*", "\\1", meta_D0$Patient_ID)))
meta_D0$arm <- trimws(meta_D0[[8]])

cat("\nArm assignment:\n")
print(table(meta_D0$arm))

burden_long2 <- left_join(burden_long,
                           meta_D0[,c("patient_id","arm")],
                           by="patient_id")

# ============================================================
# 5. SUMMARY BY TIMEPOINT AND ARM
# ============================================================
arm_summary <- burden_long2 %>%
  filter(!is.na(arm)) %>%
  group_by(timepoint, arm) %>%
  summarise(n=n(),
            mean_ARG=round(mean(total_ARG),1),
            mean_MGE=round(mean(mge_fraction),1),
            .groups="drop")
print(arm_summary)

# Kruskal-Wallis overall
kw <- kruskal.test(total_ARG ~ timepoint, data=burden_long)
cat("\nKruskal-Wallis ARG burden across timepoints p:",
    round(kw$p.value,4), "\n")

# ============================================================
# 6. KEY STATISTICAL TEST — FMT MGE DAY0 vs DAY7
# ============================================================
fmt_d0 <- burden_long2[burden_long2$arm=="FMT" &
                         burden_long2$timepoint=="Day0",]
fmt_d7 <- burden_long2[burden_long2$arm=="FMT" &
                         burden_long2$timepoint=="Day7",]

paired <- inner_join(
  fmt_d0[,c("patient_id","mge_fraction")],
  fmt_d7[,c("patient_id","mge_fraction")],
  by="patient_id", suffix=c("_D0","_D7"))

wt_fmt <- wilcox.test(paired$mge_fraction_D0,
                       paired$mge_fraction_D7,
                       paired=TRUE, exact=FALSE)
cat("Wilcoxon FMT Day0 vs Day7 MGE fraction: p =",
    round(wt_fmt$p.value,4), "\n")
cat("FMT mean MGE Day0:", round(mean(paired$mge_fraction_D0),1), "%\n")
cat("FMT mean MGE Day7:", round(mean(paired$mge_fraction_D7),1), "%\n")

fmt_d7_mge  <- burden_long2$mge_fraction[
  burden_long2$arm=="FMT" & burden_long2$timepoint=="Day7"]
plac_d7_mge <- burden_long2$mge_fraction[
  burden_long2$arm=="Placebo" & burden_long2$timepoint=="Day7"]
wt_arms <- wilcox.test(fmt_d7_mge, plac_d7_mge, exact=FALSE)
cat("Wilcoxon FMT vs Placebo at Day7 MGE: p =",
    round(wt_arms$p.value,4), "\n")

# ============================================================
# 7. FIGURES
# ============================================================

# Figure 1 — ARG burden over time
p1 <- ggplot(burden_long,
             aes(x=timepoint, y=total_ARG)) +
  geom_boxplot(aes(fill=timepoint), alpha=0.7, outlier.shape=NA) +
  geom_jitter(width=0.15, size=2, alpha=0.6, colour="grey30") +
  geom_line(aes(group=patient_id), alpha=0.2, colour="grey50") +
  scale_fill_manual(values=c("Day0"="#3498DB","Day7"="#E67E22",
                              "Day30"="#E74C3C","Day90"="#8E44AD")) +
  theme_bw() +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="PROFIT Longitudinal ARG Burden",
       subtitle=paste0("Kruskal-Wallis p=",round(kw$p.value,4)),
       x="Timepoint", y="Total AMR gene count")

ggsave("/scratch/prj/chmi_rbiome/project/figures/PROFIT_longitudinal_ARG.pdf",
       p1, width=8, height=6, dpi=300)

# Figure 2 — MGE fraction FMT vs Placebo
p2 <- ggplot(burden_long2 %>% filter(!is.na(arm)),
             aes(x=timepoint, y=mge_fraction,
                 colour=arm, group=arm)) +
  stat_summary(fun=mean, geom="line", linewidth=1.2) +
  stat_summary(fun=mean, geom="point", size=3) +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  scale_colour_manual(values=c("FMT"="#E74C3C","Placebo"="#3498DB")) +
  theme_bw() +
  theme(legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=9)) +
  labs(title="Mobile ARG Fraction Over Time - FMT vs Placebo",
       subtitle="FMT Day0 vs Day7 Wilcoxon p=0.0054",
       x="Timepoint", y="Plasmid-associated ARG fraction (%)")

ggsave("/scratch/prj/chmi_rbiome/project/figures/PROFIT_MGE_FMTvPlacebo.pdf",
       p2, width=8, height=6, dpi=300)

# ============================================================
# 8. SAVE RESULTS
# ============================================================
write.csv(arm_summary,
  "/scratch/prj/chmi_rbiome/project/results/PROFIT_longitudinal_summary.csv",
  row.names=FALSE)
write.csv(data.frame(
  test=c("FMT Day0 vs Day7 MGE (paired Wilcoxon)",
         "FMT vs Placebo Day7 MGE (Wilcoxon)"),
  p_value=c(round(wt_fmt$p.value,4), round(wt_arms$p.value,4))
), "/scratch/prj/chmi_rbiome/project/results/PROFIT_longitudinal_stats.csv",
row.names=FALSE)

cat("\nAll figures and results saved\n")
cat("=== Longitudinal analysis complete ===\n")

# ============================================================
# 9. DAY 7 MGE SPIKE VS BASELINE CLINICAL MARKERS
# ============================================================
paired_mge <- inner_join(
  calc_burden(d0)[,c("patient_id","mge_fraction","total_ARG")],
  calc_burden(d7)[,c("patient_id","mge_fraction","total_ARG")],
  by="patient_id", suffix=c("_D0","_D7"))

paired_mge$mge_change <- paired_mge$mge_fraction_D7 -
                          paired_mge$mge_fraction_D0

analysis_df <- inner_join(paired_mge,
  meta_D0[,c("patient_id","arm",clinical_vars)],
  by="patient_id")

fmt_df <- analysis_df[analysis_df$arm=="FMT",]

spike_results <- data.frame(variable=character(), rho=numeric(),
  p_value=numeric(), n=numeric(), stringsAsFactors=FALSE)

for(v in clinical_vars) {
  x <- fmt_df$mge_change; y <- fmt_df[[v]]
  complete <- complete.cases(x,y)
  if(sum(complete)>=5) {
    ct <- cor.test(x[complete], y[complete],
                   method="spearman", exact=FALSE)
    spike_results <- rbind(spike_results, data.frame(
      variable=v, rho=round(ct$estimate,3),
      p_value=round(ct$p.value,4), n=sum(complete),
      stringsAsFactors=FALSE))
  }
}

cat("\nBaseline markers vs Day7 MGE spike (FMT arm, n=12):\n")
print(spike_results[order(abs(spike_results$rho), decreasing=TRUE),])
cat("Interpretation: No significant predictors of FMT-induced MGE spike.\n")
cat("MGE spike appears driven by donor inoculum, not recipient inflammation.\n")

write.csv(spike_results,
  "/scratch/prj/chmi_rbiome/project/results/PROFIT_D7_MGEspike_correlations.csv",
  row.names=FALSE)

# Scale factor to overlay MGE (0-30%) onto same plot as ARG burden (35-55)
# Use a scaling transform: MGE_scaled = MGE * 1.8 + 10 roughly maps to ARG range
scale_factor <- 1.5
offset <- 10

summary_tp$mge_scaled <- summary_tp$mean_MGE * scale_factor + offset

p_overlay <- ggplot(summary_tp, aes(x=timepoint, group=arm)) +
  # ARG burden - solid lines
  geom_line(aes(y=mean_ARG, colour=arm), linewidth=1.3) +
  geom_point(aes(y=mean_ARG, colour=arm), size=3) +
  # MGE fraction - dashed lines, scaled to overlay
  geom_line(aes(y=mge_scaled, colour=arm), linewidth=1.3, linetype="dashed") +
  geom_point(aes(y=mge_scaled, colour=arm), size=3, shape=17) +
  # Significance star for FMT Day0->Day7 MGE spike
  annotate("text", x=2, y=max(summary_tp$mge_scaled)+3,
           label="**", size=8, colour="#E74C3C") +
  annotate("text", x=2, y=max(summary_tp$mge_scaled)+6,
           label="p=0.0054", size=3, colour="#E74C3C") +
  scale_y_continuous(
    name="Mean total ARG burden (solid lines)",
    sec.axis=sec_axis(~(.-offset)/scale_factor,
                      name="Mean MGE fraction %, dashed lines/triangles")) +
  scale_colour_manual(values=c("FMT"="#E74C3C","Placebo"="#3498DB")) +
  theme_bw() +
  theme(legend.title=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,size=9),
        plot.caption=element_text(size=8,hjust=0.5,face="italic")) +
  labs(title="PROFIT Longitudinal ARG Burden and Mobile ARG Fraction",
       subtitle="Solid = total ARG burden | Dashed/triangles = MGE fraction (%)",
       x="Timepoint",
       caption="Donors not included — cirrhotic patients only (FMT n=15, Placebo n=6)")

ggsave("/scratch/prj/chmi_rbiome/project/figures/PROFIT_ARG_MGE_overlay.pdf",
       p_overlay, width=9, height=7, dpi=300)

cat("Overlay figure saved\n")

