# ============================================================
# PROFIT Clinical Metadata Analysis
# Complete reproducible script covering:
# 1. Data import and cleaning
# 2. Normality testing
# 3. Novel correlation analysis
# 4. Visualisation
# Author: Eleni Kalopedis
# Date: 2026-05-31
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
  skip = 2, header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)

day7 <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Day7.csv",
  header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)

day30 <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Day30.csv",
  header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)

day90 <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Day90.csv",
  header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)

cat("Baseline rows:", nrow(baseline), "\n")
cat("Day7 rows:", nrow(day7), "\n")
cat("Day30 rows:", nrow(day30), "\n")
cat("Day90 rows:", nrow(day90), "\n")

# ============================================================
# 2. RENAME KEY COLUMNS — BASELINE
# ============================================================

colnames(baseline)[1]  <- "patient_id"
colnames(baseline)[2]  <- "pseudonymised_id"
colnames(baseline)[9]  <- "MELD"
colnames(baseline)[29] <- "Ammonia_plasma"
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
# 3. FORCE NUMERIC CONVERSION
# ============================================================

key_cols <- c("E_faecalis_copies", "E_coli_copies", "MELD",
              "Calprotectin", "faecal_ammonia", "Ammonia_plasma",
              "faecal_IL17A", "faecal_IL17E", "faecal_IL17F",
              "faecal_IL21", "faecal_IL22", "faecal_IFNg",
              "faecal_IL10", "faecal_IL1b", "faecal_IL6",
              "faecal_TNFa", "faecal_IL12", "faecal_IL23",
              "faecal_IL8", "faecal_Dlactate", "plasma_Dlactate")

for(col in key_cols) {
  if(col %in% colnames(baseline)) {
    baseline[[col]] <- as.numeric(as.character(baseline[[col]]))
  }
}

# ============================================================
# 4. FILTER BY TREATMENT ARM
# ============================================================

fmt     <- baseline[which(baseline$IMP == 1), ]
placebo <- baseline[which(baseline$IMP == 2), ]

cat("\nTreatment arms:\n")
cat("FMT rows:", nrow(fmt), "\n")
cat("Placebo rows:", nrow(placebo), "\n")

# ============================================================
# 5. NORMALITY TESTING — KEY VARIABLES
# ============================================================

cat("\n=== Shapiro-Wilk Normality Tests (D0) ===\n")

d0_fmt <- fmt[grep("_D0$", fmt$pseudonymised_id), ]

normality_vars <- c("faecal_ammonia", "faecal_IL22", "faecal_IL1b",
                    "MELD", "Calprotectin", "E_faecalis_copies",
                    "E_coli_copies")

for(var in normality_vars) {
  vals <- na.omit(d0_fmt[[var]])
  if(length(vals) >= 3) {
    sw <- shapiro.test(vals)
    distribution <- ifelse(sw$p.value > 0.05, "NORMAL", "NON-NORMAL")
    cat(var, "— W =", round(sw$statistic, 3),
        ", p =", round(sw$p.value, 4),
        "—", distribution, "\n")
  }
}

cat("\nConclusion: Non-parametric tests (Spearman) used throughout\n")

# ============================================================
# 6. NOVEL CORRELATION ANALYSIS
# ============================================================

cytokines <- c("faecal_IL17A", "faecal_IL17E", "faecal_IL17F",
               "faecal_IL21", "faecal_IL22", "faecal_IFNg",
               "faecal_IL10", "faecal_IL1b", "faecal_IL6",
               "faecal_TNFa", "faecal_IL12", "faecal_IL23",
               "faecal_IL8")

outcomes <- c("E_faecalis_copies", "E_coli_copies",
              "MELD", "Calprotectin", "faecal_ammonia")

results <- data.frame(
  arm = character(), timepoint = character(),
  outcome = character(), cytokine = character(),
  n = integer(), rho = numeric(),
  p_value = numeric(), stringsAsFactors = FALSE
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
                             exact = FALSE)

          results <- rbind(results, data.frame(
            arm = arm_name, timepoint = tp,
            outcome = outcome, cytokine = cytokine,
            n = nrow(cor_sub),
            rho = round(result$estimate, 3),
            p_value = round(result$p.value, 4),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

results$p_adjusted <- round(p.adjust(results$p_value,
                                     method = "BH"), 4)

# ============================================================
# 7. SAVE RESULTS
# ============================================================

write.csv(results,
  "/scratch/prj/chmi_rbiome/project/results/PROFIT_correlations_all.csv",
  row.names = FALSE)

# Print summary
cat("\n=== Correlation Analysis Summary ===\n")
cat("Total correlations:", nrow(results), "\n")
sig_fdr <- results[!is.na(results$p_adjusted) &
                   results$p_adjusted <= 0.05, ]
cat("FDR significant (p_adj <= 0.05):", nrow(sig_fdr), "\n")
exploratory <- results[!is.na(results$p_value) &
                       results$p_value <= 0.05, ]
cat("Nominally significant (p <= 0.05):", nrow(exploratory), "\n")

cat("\n=== Nominally Significant Results ===\n")
exploratory <- exploratory[order(abs(exploratory$rho),
                                 decreasing = TRUE), ]
print(exploratory)
