# ============================================================
# PROFIT Patient Metabolic Analysis
# p20 (CRKP/K. pneumoniae) vs p23 (VR-CDI/C. difficile)
# Bin-level metabolic pathway comparison
# Author: Eleni Kalopedis
# Date: 2026-06-18
# ============================================================

library(ggplot2)

fig_path <- "/scratch/prj/chmi_rbiome/project/figures"
base_met <- "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_ProxiMeta"

# ============================================================
# 1. p20 — K. pneumoniae bin_8 metabolic profile
# Key finding: bin_8 confirmed by grep Klebsiella in
# p20-d1_1180293_1094/amr_composition.tsv
# blaKPC-2 on pMAG_4 (plasmid MAG) — plasmid-borne confirmed
# ============================================================

p20_bin8 <- read.table(
  file.path(base_met,
    "p20metabolism/cluster_module_completion/bin_8.tsv"),
  header=TRUE, sep="\t", stringsAsFactors=FALSE,
  quote="", fill=TRUE)

cat("=== p20 bin_8 K. pneumoniae — complete pathways (>70%) ===\n")
p20_complete <- p20_bin8[p20_bin8$Percent_steps_found >= 70,
                          c("Module_name","Percent_steps_found")]
cat("Total complete pathways:", nrow(p20_complete), "\n")
print(p20_complete)

# Key resistance/secretion pathways
key_systems_p20 <- c(
  "Type IV secretion",
  "Type VI secretion",
  "capsule",
  "acid",
  "AcrAB",
  "MdlAB",
  "CAMP")
cat("\nKey pathways in p20 bin_8:\n")
for(s in key_systems_p20) {
  hits <- p20_bin8[grepl(s, p20_bin8$Module_name, ignore.case=TRUE),
                   c("Module_name","Percent_steps_found")]
  if(nrow(hits)>0) print(hits)
}

# ============================================================
# 2. p23 — C. difficile bin_17 metabolic profile
# Key finding: bin_17 confirmed by grep Clostridioides in
# p23-d1_1180389_1447/amr_composition.tsv
# All resistance genes on bin_17 directly (chromosomal) — NOT pMAG
# ============================================================

p23_bin17 <- read.table(
  file.path(base_met,
    "p23metabolism/p23-d1_1180389_Metabolism_Files/metabolism_results/cluster_module_completion/bin_17.tsv"),
  header=TRUE, sep="\t", stringsAsFactors=FALSE,
  quote="", fill=TRUE)

cat("\n=== p23 bin_17 C. difficile — complete pathways (>70%) ===\n")
p23_complete <- p23_bin17[p23_bin17$Percent_steps_found >= 70,
                           c("Module_name","Percent_steps_found")]
cat("Total complete pathways:", nrow(p23_complete), "\n")
print(p23_complete)

# ============================================================
# 3. COMPARATIVE FIGURE — Mechanistic basis of differential
#    FMT efficacy
# ============================================================
comparison <- data.frame(
  pathway=c("Type IV secretion\n(plasmid transfer)",
            "Type VI secretion\n(competitor killing)",
            "Capsule synthesis\n(immune evasion)",
            "Acid tolerance\n(EvgS-EvgA)",
            "AcrAB-TolC efflux\n(multidrug resistance)",
            "Sporulation\n(persistence)",
            "Vancomycin resistance\nregulatory (VanB/VanE)"),
  Kpneumoniae=c(83.3,100,100,100,100,0,0),
  Cdifficile=c(0,0,0,0,50,14.3,100)
)

df_long <- rbind(
  data.frame(pathway=comparison$pathway,
             organism="K. pneumoniae bin_8\n(p20 CRKP)",
             completeness=comparison$Kpneumoniae),
  data.frame(pathway=comparison$pathway,
             organism="C. difficile bin_17\n(p23 VR-CDI)",
             completeness=comparison$Cdifficile)
)

p <- ggplot(df_long,
            aes(x=reorder(pathway,-completeness),
                y=completeness, fill=organism)) +
  geom_bar(stat="identity", position="dodge", alpha=0.85) +
  geom_hline(yintercept=70, linetype="dashed",
             colour="grey40", linewidth=0.5) +
  annotate("text", x=6.5, y=72,
           label="70% threshold", size=3, colour="grey40") +
  scale_fill_manual(values=c(
    "K. pneumoniae bin_8\n(p20 CRKP)"="#E74C3C",
    "C. difficile bin_17\n(p23 VR-CDI)"="#3498DB")) +
  coord_flip() +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position="bottom",
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,size=8)) +
  labs(title="Metabolic Basis of Differential FMT Efficacy",
       subtitle=paste0("K. pneumoniae retains HGT and persistence capacity;",
                       " C. difficile does not"),
       x="", y="Pathway completeness (%)")

ggsave(file.path(fig_path,"PROFIT_p20_vs_p23_metabolic.pdf"),
       p, width=10, height=7, dpi=300)
cat("Figure saved\n")

