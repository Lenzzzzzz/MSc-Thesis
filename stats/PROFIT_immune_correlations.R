# ============================================================
# PROFIT Immune Marker × ARG Correlations
# Faecal cytokines, plasma cytokines, barrier markers,
# immunoglobulins, endotoxin, bacterial copy numbers
# Author: Eleni Kalopedis
# Date: 2026-06-17
# ============================================================

library(dplyr)
library(ggplot2)

base     <- "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_ProxiMeta"
fig_path <- "/scratch/prj/chmi_rbiome/project/figures"
res_path <- "/scratch/prj/chmi_rbiome/project/results"

# ============================================================
# 1. LOAD ARG BURDEN
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
# ============================================================
immune_raw <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Immune_Metadata.csv",
  stringsAsFactors=FALSE, check.names=FALSE, header=FALSE)
colnames(immune_raw) <- paste0("V",1:ncol(immune_raw))
immune_values <- immune_raw[4:nrow(immune_raw),]

# Faecal cytokines (normalised columns, odd positions 41-65)
cytokine_names <- c("IL17A","IL17E","IL17F","IL21","IL22",
                    "IFNg","IL10","IL1b","IL6","TNFa",
                    "IL12","IL23","IL8","FABP2_faecal")
cytokine_cols  <- seq(41,67,by=2)

# Barrier markers
# FABP2 faecal=67, plasma=69
# Dlactate faecal=71, plasma=73

# Extended markers
# Endotoxin=77, Efaecalis=79, Ecoli=80
# IgA faecal=84, IgG faecal=88, IgM faecal=92
# Plasma cytokines: IL17E=99,IL17F=101,IL21=103,IL22=105,
#                   IFNg=107,IL10=109,IL1b=111,IL6=113,TNFa=115,IL12=117

immune_df <- data.frame(
  Age=suppressWarnings(as.numeric(immune_values[["V2"]])),
  MELD=suppressWarnings(as.numeric(immune_values[["V6"]])),
  FABP2_plasma=suppressWarnings(as.numeric(immune_values[["V69"]])),
  Dlactate_faecal=suppressWarnings(as.numeric(immune_values[["V71"]])),
  Dlactate_plasma=suppressWarnings(as.numeric(immune_values[["V73"]])),
  Endotoxin=suppressWarnings(as.numeric(immune_values[["V77"]])),
  Efaecalis_copies=suppressWarnings(as.numeric(immune_values[["V79"]])),
  Ecoli_copies=suppressWarnings(as.numeric(immune_values[["V80"]])),
  IgA_faecal=suppressWarnings(as.numeric(immune_values[["V84"]])),
  IgG_faecal=suppressWarnings(as.numeric(immune_values[["V88"]])),
  IgM_faecal=suppressWarnings(as.numeric(immune_values[["V92"]])),
  pIL17E=suppressWarnings(as.numeric(immune_values[["V99"]])),
  pIL17F=suppressWarnings(as.numeric(immune_values[["V101"]])),
  pIL21=suppressWarnings(as.numeric(immune_values[["V103"]])),
  pIL22=suppressWarnings(as.numeric(immune_values[["V105"]])),
  pIFNg=suppressWarnings(as.numeric(immune_values[["V107"]])),
  pIL10=suppressWarnings(as.numeric(immune_values[["V109"]])),
  pIL1b=suppressWarnings(as.numeric(immune_values[["V111"]])),
  pIL6=suppressWarnings(as.numeric(immune_values[["V113"]])),
  pTNFa=suppressWarnings(as.numeric(immune_values[["V115"]])),
  pIL12=suppressWarnings(as.numeric(immune_values[["V117"]])),
  stringsAsFactors=FALSE)

for(i in seq_along(cytokine_cols)) {
  col_name <- paste0("V",cytokine_cols[i])
  immune_df[[cytokine_names[i]]] <- suppressWarnings(
    as.numeric(immune_values[[col_name]]))
}

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

immune_df$MELD_num <- immune_df$MELD
immune_df$Age_num  <- immune_df$Age
joined <- merge(immune_df,
                meta_D0[,c("patient_id","MELD_num","Age_num")],
                by=c("MELD_num","Age_num"))
joined <- joined[!duplicated(joined$patient_id),]
final_df <- merge(joined, burden, by="patient_id")
cat("Final dataset:", nrow(final_df), "patients\n")

# ============================================================
# 4. SPEARMAN CORRELATIONS
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

all_immune_vars <- c(cytokine_names,
                     "FABP2_plasma","Dlactate_faecal","Dlactate_plasma",
                     "Endotoxin","Efaecalis_copies","Ecoli_copies",
                     "IgA_faecal","IgG_faecal","IgM_faecal",
                     "pIL17E","pIL17F","pIL21","pIL22","pIFNg",
                     "pIL10","pIL1b","pIL6","pTNFa","pIL12")

arg_results  <- run_cors(final_df, "total_ARG",    all_immune_vars)
mge_results  <- run_cors(final_df, "mge_fraction", all_immune_vars)

cat("\n=== Significant immune markers vs ARG burden (p<0.1) ===\n")
print(arg_results[arg_results$p_value<0.1,])
cat("\n=== Significant immune markers vs MGE fraction (p<0.1) ===\n")
print(mge_results[mge_results$p_value<0.1,])

# ============================================================
# 5. SUMMARY HEATMAP — FIGURE C.2
# ============================================================
key_vars <- c("GGT","Ammonia","Calprotectin",
              "Endotoxin","FABP2_faecal","IgM_faecal",
              "IL12","IL23","IL8","pIL22","pIL10","IL10","pIL1b")

# Use previously computed values
heatmap_df <- data.frame(
  variable=rep(key_vars,2),
  outcome=c(rep("Total ARG burden",length(key_vars)),
             rep("MGE fraction",length(key_vars))),
  rho=c(-0.765,-0.755,-0.426,-0.078,-0.078,-0.623,
        -0.551,0.517,0.488,0.461,
        0.561,0.530,0.523,
        -0.277,-0.235,0.294,-0.532,0.441,-0.532,
        -0.352,-0.184,0.123,0.152,
        0.090,0.135,0.189),
  p_value=c(0.0003,0.0005,0.088,0.765,0.765,0.0076,
            0.022,0.034,0.047,0.063,
            0.019,0.029,0.055,
            0.282,0.363,0.252,0.028,0.076,0.028,
            0.166,0.480,0.639,0.560,
            0.733,0.606,0.517),
  stringsAsFactors=FALSE)

heatmap_df$sig <- ifelse(heatmap_df$p_value<0.05,"*",
                   ifelse(heatmap_df$p_value<0.1,"+",""))
heatmap_df$variable <- factor(heatmap_df$variable,
  levels=rev(c("GGT","Ammonia","Calprotectin",
               "Endotoxin","FABP2_faecal","IgM_faecal",
               "IL12","IL23","IL8","pIL22","pIL10","IL10","pIL1b")))

p_heatmap <- ggplot(heatmap_df,
  aes(x=outcome, y=variable, fill=rho)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=paste0(sprintf("%.2f",rho),"\n",sig)),
            size=3.2) +
  scale_fill_gradient2(low="#3498DB",mid="white",high="#E74C3C",
                        midpoint=0,limits=c(-0.9,0.9),
                        name="Spearman rho") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,size=8),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=9)) +
  labs(title="Immune-AMR Axis: Summary Spearman Correlations",
       subtitle="* p<0.05   + p<0.10",
       x="",y="Immune/Clinical Marker")

ggsave(file.path(fig_path,"PROFIT_immune_ARG_summary_heatmap.pdf"),
       p_heatmap, width=9, height=8, dpi=300)

# ============================================================
# 6. SAVE ALL RESULTS
# ============================================================
write.csv(arg_results,
  file.path(res_path,"PROFIT_complete_immune_ARG_correlations.csv"),
  row.names=FALSE)
write.csv(mge_results,
  file.path(res_path,"PROFIT_complete_immune_MGE_correlations.csv"),
  row.names=FALSE)
cat("All complete immune correlations saved\n")

