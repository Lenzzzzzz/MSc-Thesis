# ============================================================
# PROFIT Correlation Analysis
# Novel correlations not covered in pre-print
# Author: Eleni Kalopedis
# Date: 2026-05-30
# ============================================================

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# ============================================================
# 1. IMPORT DATA
# ============================================================

baseline <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Baseline.csv.csv",
  skip = 2,
  header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)

# ============================================================
# 2. RENAME KEY COLUMNS
# ============================================================

colnames(baseline)[1]  <- "patient_id"
colnames(baseline)[2]  <- "pseudonymised_id"
colnames(baseline)[9]  <- "MELD"
colnames(baseline)[29] <- "Ammonia"
colnames(baseline)[42] <- "Calprotectin"
colnames(baseline)[44] <- "faecal_IL17A"
colnames(baseline)[46] <- "faecal_IL17E"
colnames(baseline)[48] <- "faecal_IL17F"
colnames(baseline)[50] <- "faecal_IL21"
colnames(baseline)[52] <- "faecal_IL22"
colnames(baseline)[54] <- "faecal_IFNg"
colnames(baseline)[56] <- "faecal_IL10"
colnames(baseline)[58] <- "faecal_IL1b"
colnames(baseline)[60] <- "faecal_IL6"
colnames(baseline)[62] <- "faecal_TNFa"
colnames(baseline)[64] <- "faecal_IL12"
colnames(baseline)[66] <- "faecal_IL23"
colnames(baseline)[68] <- "faecal_IL8"
colnames(baseline)[74] <- "faecal_Dlactate"
colnames(baseline)[76] <- "plasma_Dlactate"
colnames(baseline)[78] <- "faecal_ammonia"
colnames(baseline)[82] <- "E_faecalis_copies"
colnames(baseline)[83] <- "E_coli_copies"

# ============================================================
# 3. FILTER BY TREATMENT ARM
# ============================================================
# Force outcome columns to numeric
for(col in c("E_faecalis_copies", "E_coli_copies",
             "MELD", "Calprotectin", "faecal_ammonia",
             "faecal_IL17A", "faecal_IL17E", "faecal_IL17F",
             "faecal_IL21", "faecal_IL22", "faecal_IFNg",
             "faecal_IL10", "faecal_IL1b", "faecal_IL6",
             "faecal_TNFa", "faecal_IL12", "faecal_IL23",
             "faecal_IL8")) {
  baseline[[col]] <- as.numeric(as.character(baseline[[col]]))
}
fmt     <- baseline[which(baseline$IMP == 1), ]
placebo <- baseline[which(baseline$IMP == 2), ]

# ============================================================
# 4. DEFINE VARIABLES FOR CORRELATION
# ============================================================

cytokines <- c("faecal_IL17A", "faecal_IL17E", "faecal_IL17F",
               "faecal_IL21", "faecal_IL22", "faecal_IFNg",
               "faecal_IL10", "faecal_IL1b", "faecal_IL6",
               "faecal_TNFa", "faecal_IL12", "faecal_IL23",
               "faecal_IL8")

outcomes <- c("E_faecalis_copies", "E_coli_copies",
              "MELD", "Calprotectin", "faecal_ammonia")

# ============================================================
# 5. RUN CORRELATIONS
# ============================================================

results <- data.frame(
  arm        = character(),
  timepoint  = character(),
  outcome    = character(),
  cytokine   = character(),
  n          = integer(),
  rho        = numeric(),
  p_value    = numeric(),
  stringsAsFactors = FALSE
)

for(arm_name in c("FMT", "Placebo")) {
  arm_data <- if(arm_name == "FMT") fmt else placebo

  for(tp in c("D0", "D7", "D30", "D90")) {
    tp_data <- arm_data[grep(paste0("_", tp, "$"),
                             arm_data$pseudonymised_id), ]

    for(outcome in outcomes) {
      for(cytokine in cytokines) {

        if(!outcome %in% colnames(tp_data) |
           !cytokine %in% colnames(tp_data)) next

        cor_sub <- na.omit(tp_data[, c(outcome, cytokine)])

        if(nrow(cor_sub) > 4) {
          result <- cor.test(cor_sub[[outcome]],
                             cor_sub[[cytokine]],
                             method = "spearman",
                             exact  = FALSE)

          results <- rbind(results, data.frame(
            arm       = arm_name,
            timepoint = tp,
            outcome   = outcome,
            cytokine  = cytokine,
            n         = nrow(cor_sub),
            rho       = round(result$estimate, 3),
            p_value   = round(result$p.value,  4),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

# ============================================================
# 6. APPLY FDR CORRECTION
# ============================================================

results$p_adjusted <- round(p.adjust(results$p_value, method = "BH"), 4)

# ============================================================
# 7. SAVE RESULTS
# ============================================================

write.csv(results,
  "/scratch/prj/chmi_rbiome/project/results/PROFIT_correlations_all.csv",
  row.names = FALSE)

# Print significant results
sig <- results[results$p_adjusted <= 0.05, ]
cat("\n=== Significant correlations (FDR <= 0.05) ===\n")
cat("Total:", nrow(sig), "\n\n")
print(sig)

# Print top 10 by rho magnitude regardless of significance
cat("\n=== Top 10 correlations by rho magnitude ===\n")
results_sorted <- results[order(abs(results$rho), decreasing = TRUE), ]
print(head(results_sorted, 10))

