# ============================================================
# PROFIT Immune Marker x ARG Correlations
# Faecal/plasma cytokines, barrier markers, immunoglobulins,
# endotoxin, clinical markers (GGT)
# CORRECTED 2026-06-20: plasma cytokine columns were
# previously misaligned (off by 1-2 positions from IL-21
# onward) and the summary heatmap previously used HARDCODED
# values disconnected from the live computation. This version
# fixes both issues - all values now come directly from
# run_cors() with verified column mappings.
# Author: Eleni Kalopedis
# ============================================================

library(dplyr)
library(ggplot2)

base     <- "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_ProxiMeta"
fig_path <- "/scratch/prj/chmi_rbiome/project/figures"
res_path <- "/scratch/prj/chmi_rbiome/project/results"

# ============================================================
# 1. LOAD ARG BURDEN (DAY 0, DONORS EXCLUDED)
# ============================================================
import_timepoint <- function(timepoint_path, timepoint_label) {
  folders <- list.dirs(timepoint_path, recursive=FALSE)
  amr_list <- lapply(folders, function(f) {
    amr_file <- file.path(f, "amr_composition.tsv")
    if(file.exists(amr_file)) {
      df <- tryCatch(read.table(amr_file, header=TRUE, sep="\t",
                     stringsAsFactors=FALSE, quote=""),
                     error=function(e) NULL)
      if(!is.null(df)) {
        df$folder <- basename(f)
        df$timepoint <- timepoint_label
        df$patient_id <- gsub("(p\\d+)-.*","\\1", basename(f))
        df
      }
    }
  })
  do.call(rbind, amr_list[!sapply(amr_list, is.null)])
}

d0 <- import_timepoint(file.path(base,"PME_PROFIT"),"Day0")
d0 <- d0[!grepl("^d", d0$folder),]
amr_d0 <- d0[d0$Element.type=="AMR",]
burden <- amr_d0 %>%
  group_by(patient_id) %>%
  summarise(total_ARG=n(),
            plasmid_ARGs=sum(Gene.contig.type=="plasmid",na.rm=TRUE),
            mge_fraction=round(plasmid_ARGs/total_ARG*100,1),
            .groups="drop")

# ============================================================
# 2. LOAD AND PARSE IMMUNE METADATA
# VERIFIED COLUMN MAP (text-match confirmed 2026-06-20):
# Faecal cytokines: IL17A=41,IL17E=43,IL17F=45,IL21=47,IL22=49,
#   IFNg=51,IL10=53,IL1b=55,IL6=57,TNFa=59,IL12=61,IL23=63,IL8=65
# Barrier: FABP2 faecal=67,plasma=69; Dlactate faecal=71,plasma=73
# Analytes: Ammonia=75, Endotoxin=77(F)/125(P)
# Bacterial copies: Efaecalis=79, Ecoli=80
# Immunoglobulins: IgA faecal=84/plasma=86, IgG faecal=88/plasma=90,
#   IgM faecal=92/plasma=94
# Plasma cytokines (CORRECTED, irregular spacing):
#   IL17A=96,IL17E=99,IL17F=102,IL21=105,IL22=107,IFNg=109,
#   IL10=111,IL1b=113,IL6=115,TNFa=117,IL12=119,IL23=121,IL8=123
# ============================================================
immune_raw <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Immune_Metadata.csv",
  stringsAsFactors=FALSE, check.names=FALSE, header=FALSE)
colnames(immune_raw) <- paste0("V",1:ncol(immune_raw))
immune_values <- immune_raw[4:nrow(immune_raw),]

immune_df <- data.frame(
  MELD=suppressWarnings(as.numeric(immune_values[["V6"]])),
  Age=suppressWarnings(as.numeric(immune_values[["V2"]])),
  IL12_faecal=suppressWarnings(as.numeric(immune_values[["V61"]])),
  IL23_faecal=suppressWarnings(as.numeric(immune_values[["V63"]])),
  IL8_faecal=suppressWarnings(as.numeric(immune_values[["V65"]])),
  IL10_faecal=suppressWarnings(as.numeric(immune_values[["V53"]])),
  IL22_faecal=suppressWarnings(as.numeric(immune_values[["V49"]])),
  IFNg_faecal=suppressWarnings(as.numeric(immune_values[["V51"]])),
  FABP2_faecal=suppressWarnings(as.numeric(immune_values[["V67"]])),
  FABP2_plasma=suppressWarnings(as.numeric(immune_values[["V69"]])),
  Dlactate_faecal=suppressWarnings(as.numeric(immune_values[["V71"]])),
  Dlactate_plasma=suppressWarnings(as.numeric(immune_values[["V73"]])),
  Ammonia_faecal=suppressWarnings(as.numeric(immune_values[["V75"]])),
  Endotoxin_faecal=suppressWarnings(as.numeric(immune_values[["V77"]])),
  Endotoxin_plasma=suppressWarnings(as.numeric(immune_values[["V125"]])),
  Efaecalis_copies=suppressWarnings(as.numeric(immune_values[["V79"]])),
  Ecoli_copies=suppressWarnings(as.numeric(immune_values[["V80"]])),
  IgA_faecal=suppressWarnings(as.numeric(immune_values[["V84"]])),
  IgA_plasma=suppressWarnings(as.numeric(immune_values[["V86"]])),
  IgG_faecal=suppressWarnings(as.numeric(immune_values[["V88"]])),
  IgG_plasma=suppressWarnings(as.numeric(immune_values[["V90"]])),
  IgM_faecal=suppressWarnings(as.numeric(immune_values[["V92"]])),
  IgM_plasma=suppressWarnings(as.numeric(immune_values[["V94"]])),
  pIL17A=suppressWarnings(as.numeric(immune_values[["V96"]])),
  pIL17E=suppressWarnings(as.numeric(immune_values[["V99"]])),
  pIL17F=suppressWarnings(as.numeric(immune_values[["V102"]])),
  pIL21=suppressWarnings(as.numeric(immune_values[["V105"]])),
  pIL22=suppressWarnings(as.numeric(immune_values[["V107"]])),
  pIFNg=suppressWarnings(as.numeric(immune_values[["V109"]])),
  pIL10=suppressWarnings(as.numeric(immune_values[["V111"]])),
  pIL1b=suppressWarnings(as.numeric(immune_values[["V113"]])),
  pIL6=suppressWarnings(as.numeric(immune_values[["V115"]])),
  pTNFa=suppressWarnings(as.numeric(immune_values[["V117"]])),
  pIL12=suppressWarnings(as.numeric(immune_values[["V119"]])),
  pIL23=suppressWarnings(as.numeric(immune_values[["V121"]])),
  pIL8=suppressWarnings(as.numeric(immune_values[["V123"]])),
  stringsAsFactors=FALSE)

# ============================================================
# 3. JOIN WITH PATIENT IDs
# ============================================================
meta_raw2 <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Baseline.csv",
  stringsAsFactors=FALSE, check.names=FALSE,
  header=FALSE, skip=2)
colnames(meta_raw2) <- as.character(meta_raw2[1,])
meta_raw2 <- meta_raw2[-1,]
colnames(meta_raw2)[1] <- "Patient_ID"
colnames(meta_raw2) <- make.unique(colnames(meta_raw2), sep="_")
meta_D0 <- meta_raw2[grepl("D0$", meta_raw2$Patient_ID),]
meta_D0$patient_id <- paste0("p", as.integer(
  gsub("P0100(\\d+).*","\\1", meta_D0$Patient_ID)))
meta_D0$MELD_num <- suppressWarnings(as.numeric(meta_D0$MELD))
meta_D0$Age_num  <- suppressWarnings(as.numeric(
  meta_D0$`Age at baseline/years`))
meta_D0$GGT_num <- suppressWarnings(as.numeric(meta_D0$GGT))

immune_df$MELD_num <- immune_df$MELD
immune_df$Age_num  <- immune_df$Age
joined <- merge(immune_df,
                meta_D0[,c("patient_id","MELD_num","Age_num","GGT_num")],
                by=c("MELD_num","Age_num"))
joined <- joined[!duplicated(joined$patient_id),]
final_df <- merge(joined, burden, by="patient_id")
cat("Final dataset:", nrow(final_df), "patients\n")

# ============================================================
# 4. SPEARMAN CORRELATIONS - ALL LIVE COMPUTATION
# ============================================================
run_cors <- function(df, y_var, predictors) {
  results <- data.frame(variable=character(), rho=numeric(),
    p_value=numeric(), n=numeric(), stringsAsFactors=FALSE)
  for(v in predictors) {
    x <- df[[y_var]]; y <- df[[v]]
    complete <- complete.cases(x,y)
    if(sum(complete)>=5) {
      ct <- cor.test(x[complete],y[complete],
                     method="spearman",exact=FALSE)
      results <- rbind(results, data.frame(
        variable=v, rho=round(ct$estimate,3),
        p_value=round(ct$p.value,4), n=sum(complete),
        stringsAsFactors=FALSE))
    }
  }
  results[order(abs(results$rho),decreasing=TRUE),]
}

all_immune_vars <- c("IL12_faecal","IL23_faecal","IL8_faecal","IL10_faecal",
                     "IL22_faecal","IFNg_faecal","FABP2_faecal","FABP2_plasma",
                     "Dlactate_faecal","Dlactate_plasma","Ammonia_faecal",
                     "Endotoxin_faecal","Endotoxin_plasma",
                     "Efaecalis_copies","Ecoli_copies",
                     "IgA_faecal","IgA_plasma","IgG_faecal","IgG_plasma",
                     "IgM_faecal","IgM_plasma",
                     "pIL17A","pIL17E","pIL17F","pIL21","pIL22","pIFNg",
                     "pIL10","pIL1b","pIL6","pTNFa","pIL12","pIL23","pIL8",
                     "GGT_num")

arg_results <- run_cors(final_df, "total_ARG", all_immune_vars)
mge_results <- run_cors(final_df, "mge_fraction", all_immune_vars)

cat("\n=== ALL markers vs ARG burden, sorted by p-value ===\n")
print(arg_results[order(arg_results$p_value),])
cat("\n=== ALL markers vs MGE fraction, sorted by p-value ===\n")
print(mge_results[order(mge_results$p_value),])

write.csv(arg_results,
  file.path(res_path,"PROFIT_complete_immune_ARG_correlations_VERIFIED.csv"),
  row.names=FALSE)
write.csv(mge_results,
  file.path(res_path,"PROFIT_complete_immune_MGE_correlations_VERIFIED.csv"),
  row.names=FALSE)

# ============================================================
# 5. SUMMARY HEATMAP - 3 PANELS, LIVE VALUES ONLY (NO HARDCODING)
# Postdoc requested: split Clinical/Cytokines/Barrier panels,
# add units, add F/P labels, clarify correlation framing,
# add D-lactate. Supervisor requested: integrate plasma data.
# ============================================================
heatmap_vars <- c("GGT_num","Ammonia_faecal","IL12_faecal","IL23_faecal",
                   "IL8_faecal","Endotoxin_faecal","FABP2_faecal",
                   "Dlactate_faecal","Dlactate_plasma","IgM_faecal",
                   "IgM_plasma","pIL21","pIFNg","pIL10")

arg_sub <- arg_results[arg_results$variable %in% heatmap_vars,]
arg_sub$outcome <- "Total ARG burden"
mge_sub <- mge_results[mge_results$variable %in% heatmap_vars,]
mge_sub$outcome <- "MGE fraction"
heatmap_df <- rbind(arg_sub, mge_sub)
heatmap_df$sig <- ifelse(heatmap_df$p_value<0.05,"*",
                   ifelse(heatmap_df$p_value<0.1,"+",""))

heatmap_df$panel <- case_when(
  heatmap_df$variable %in% c("GGT_num","Ammonia_faecal") ~ "Clinical",
  heatmap_df$variable %in% c("IL12_faecal","IL23_faecal","IL8_faecal",
                              "pIL21","pIFNg","pIL10") ~ "Cytokines",
  TRUE ~ "Barrier/Mucosal")

label_map_units <- c(
  "GGT_num"="GGT (Clinical, IU/L)",
  "Ammonia_faecal"="Ammonia (F, unit TBC)",
  "IL12_faecal"="IL-12 (F, pg/mg protein*)",
  "IL23_faecal"="IL-23 (F, pg/mg protein*)",
  "IL8_faecal"="IL-8 (F, pg/mg protein*)",
  "pIL21"="IL-21 (P, pg/mL*)",
  "pIFNg"="IFNg (P, pg/mL*)",
  "pIL10"="IL-10 (P, pg/mL*)",
  "Endotoxin_faecal"="Endotoxin (F, unit TBC)",
  "FABP2_faecal"="FABP2 (F, ng/mg protein)",
  "Dlactate_faecal"="D-lactate (F, nMOL)",
  "Dlactate_plasma"="D-lactate (P, nMOL)",
  "IgM_faecal"="IgM (F, unit TBC)",
  "IgM_plasma"="IgM (P, unit TBC)"
)
heatmap_df$label <- label_map_units[heatmap_df$variable]
heatmap_df$panel <- factor(heatmap_df$panel,
                            levels=c("Clinical","Cytokines","Barrier/Mucosal"))

p_heatmap_final <- ggplot(heatmap_df,
  aes(x=outcome, y=reorder(label,rho), fill=rho)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=paste0(sprintf("%.2f",rho),sig)), size=3) +
  facet_grid(panel~., scales="free_y", space="free_y") +
  scale_fill_gradient2(low="#3498DB", mid="white", high="#E74C3C",
                        midpoint=0, limits=c(-0.8,0.8),
                        name="Spearman\nrho") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=8),
        plot.caption=element_text(size=7, hjust=0, face="italic"),
        strip.text.y=element_text(face="bold", angle=0),
        axis.text.x=element_text(size=9)) +
  labs(title="Immune, Clinical and Barrier Marker Correlations with the Resistome",
       subtitle="Each marker correlated independently against Total ARG burden and MGE fraction. * p<0.05, + p<0.1",
       x="", y="",
       caption="(F) = Faecal, (P) = Plasma. *Unit inferred, not independently confirmed.")

ggsave(file.path(fig_path,"PROFIT_immune_ARG_summary_heatmap_FINAL.pdf"),
       p_heatmap_final, width=9, height=9, dpi=300)

write.csv(heatmap_df,
  file.path(res_path,"PROFIT_complete_immune_correlations_VERIFIED.csv"),
  row.names=FALSE)

cat("All complete - script fully corrected, no hardcoded values\n")
