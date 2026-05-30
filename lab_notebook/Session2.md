# Computational Lab Notebook
## MSc Thesis â€” Investigating the 'Immune-AMR' Axis
**Author:** Eleni Andria Kalopedis
**Supervisor:** Dr. Lindsey Ann Edwards, King's College London
**Session:** 2
**Date:** 30 May 2026
**Platform:** CREATE HPC (King's College London) + Chromebook/Linux local

---

## Session Overview

Continuation of Session 1. Primary objectives were to write and submit the Trimmomatic pipeline script, set up tmux for persistent sessions, transfer and analyse the PROFIT clinical metadata in R, and run the first novel correlation analyses not covered in the pre-print.

---

## 1. HPC Portal Outage and Resolution

### What happened
The CREATE MFA portal (portal.er.kcl.ac.uk) went down mid-session returning a 50x server error for approximately 2 hours. This prevented re-login to the HPC after terminal disconnections.

### Resolution
Portal restored by KCL eResearch team. No action required on user end.

### Key learning
Always work inside a tmux session to prevent losing work during portal outages or laptop disconnections. Set up in this session â€” see Section 5.

---

## 2. Trimmomatic Script â€” First Submission

### Script written
File: `/scratch/prj/chmi_rbiome/project/scripts/MSc-Thesis/pipeline/trimmomatic.sh`

### What the script does
- Loads Trimmomatic module v0.39
- Loops through all 10 POIROT samples automatically
- Runs paired-end trimming on R1 and R2 files for each sample
- Removes TruSeq3 adapters, low quality bases (LEADING:3, TRAILING:3), and reads under 36bp (MINLEN:36)
- Saves four output files per sample (paired R1, unpaired R1, paired R2, unpaired R2)
- Logs progress to `/scratch/prj/chmi_rbiome/project/logs/`

### SLURM parameters (first submission)
```
--job-name=trimmomatic
--cpus-per-task=4
--mem=16G
--time=04:00:00
```

### Job ID: 34324675

### Issue encountered
At 50 minutes runtime, only 3 of 10 samples had completed. Estimated total runtime ~5 hours â€” exceeded the 4-hour time limit. Job cancelled manually:
```bash
scancel 34324675
```

### Output before cancellation
Three samples completed successfully:
- `KINGCO-007-AKQ014A` â€” 4 output files saved âœ…
- `KINGCO-007-AKQ014B` â€” 4 output files saved âœ…
- `KINGCO-007-AKQ01A` â€” 4 output files saved âœ…

All output in: `/scratch/prj/chmi_rbiome/project/results/trimmomatic/`

---

## 3. Trimmomatic Script â€” Resubmission

### Changes made to script
1. Time limit increased from `04:00:00` to `12:00:00`
2. Skip condition added to avoid reprocessing completed samples:

```bash
if [ -f "$OUT/${SAMPLE}_R1_paired.fastq.gz" ]; then
    echo "Skipping $SAMPLE - already done"
    continue
fi
```

### Job ID: 34326136
### Status at end of session: Running â€” 11 minutes in, 7 remaining samples

### Decision made
12-hour time limit set as conservative upper bound. At ~30 minutes per sample, 7 remaining samples = ~3.5 hours. Well within new limit.

---

## 4. PROFIT Clinical Metadata â€” Transfer to HPC

### Files transferred
Local path on Chromebook: `/mnt/chromeos/MyFiles/Downloads/LE Lab/`

Note: Files were not in the standard Linux Downloads path â€” required sharing the "LE Lab" folder with Linux via Chrome OS Files app first.

```bash
scp "/mnt/chromeos/MyFiles/Downloads/LE Lab/PROFIT_Baseline.csv.csv" k25118483@hpc.create.kcl.ac.uk:/scratch/prj/chmi_rbiome/project/metadata/
scp "/mnt/chromeos/MyFiles/Downloads/LE Lab/PROFIT_Day7.csv" k25118483@hpc.create.kcl.ac.uk:/scratch/prj/chmi_rbiome/project/metadata/
scp "/mnt/chromeos/MyFiles/Downloads/LE Lab/PROFIT_Day30.csv" k25118483@hpc.create.kcl.ac.uk:/scratch/prj/chmi_rbiome/project/metadata/
scp "/mnt/chromeos/MyFiles/Downloads/LE Lab/PROFIT_Day90.csv" k25118483@hpc.create.kcl.ac.uk:/scratch/prj/chmi_rbiome/project/metadata/
```

Note: Baseline file had double extension `PROFIT_Baseline.csv.csv` â€” carried through from Google Sheets export. No impact on functionality.

### Confirmed on HPC
```bash
ls /scratch/prj/chmi_rbiome/project/metadata/
# Output: PROFIT_Baseline.csv.csv  PROFIT_Day30.csv  PROFIT_Day7.csv  PROFIT_Day90.csv
```

---

## 5. tmux Setup

### What tmux is
A persistent terminal multiplexer â€” sessions survive laptop disconnections and screen closure.

### Setup commands
```bash
tmux new -s thesis
```

### How to reconnect after disconnection
```bash
tmux attach -t thesis
```

### Key shortcuts
| Shortcut | Action |
|---|---|
| Ctrl + B then D | Detach (leave running in background) |
| Ctrl + B then C | New window |
| Ctrl + B then 0/1/2 | Switch windows |

---

## 6. PROFIT Correlation Analysis â€” R

### Import and cleaning

```r
baseline <- read.csv(
  "/scratch/prj/chmi_rbiome/project/metadata/PROFIT_Baseline.csv.csv",
  skip = 2,
  header = TRUE,
  na.strings = c("NA", "", "#DIV/0!", "NSA", "ni")
)
```

**Dataset dimensions:** 98 rows Ã— 153 columns

**Key issue:** Original Excel file had merged headers spanning multiple columns â€” R read these as `X`, `X.1`, `X.2` etc. Resolved by manually renaming key columns by position.

### Columns renamed

| Column number | New name |
|---|---|
| 1 | patient_id |
| 2 | pseudonymised_id |
| 9 | MELD |
| 29 | Ammonia (plasma) |
| 42 | Calprotectin |
| 44 | faecal_IL17A |
| 46 | faecal_IL17E |
| 48 | faecal_IL17F |
| 50 | faecal_IL21 |
| 52 | faecal_IL22 |
| 54 | faecal_IFNg |
| 56 | faecal_IL10 |
| 58 | faecal_IL1b |
| 60 | faecal_IL6 |
| 62 | faecal_TNFa |
| 64 | faecal_IL12 |
| 66 | faecal_IL23 |
| 68 | faecal_IL8 |
| 74 | faecal_Dlactate |
| 76 | plasma_Dlactate |
| 78 | faecal_ammonia |
| 82 | E_faecalis_copies |
| 83 | E_coli_copies |

### Treatment arms confirmed
- IMP = 1: FMT group â€” 60 rows (15 patients Ã— 4 timepoints) âœ…
- IMP = 2: Placebo group â€” 24 rows (6 patients Ã— 4 timepoints) âœ…

---

## 7. Priority Correlation â€” Faecal Ammonia vs Faecal IL-22

### Rationale
Joseph's email and the pre-print both identified ammonia and IL-22 as priority markers. Initial test used plasma ammonia (column 29) â€” no significant correlation found. Switched to faecal ammonia (column 78) which is the biologically appropriate compartment for a mucosal cytokine.

### Normality testing
Both faecal ammonia and faecal IL-22 were strongly non-normal (Shapiro-Wilk p < 0.05 at all timepoints). **Spearman's rank correlation used throughout.**

### Results by timepoint (all patients combined)

| Timepoint | n | rho | p-value |
|---|---|---|---|
| D0 | 21 | +0.171 | 0.457 |
| D7 | 20 | +0.024 | 0.920 |
| D30 | 21 | +0.109 | 0.637 |
| D90 | 20 | -0.335 | 0.148 |

### Results by arm and timepoint

**FMT arm:**

| Timepoint | n | rho | p-value |
|---|---|---|---|
| D0 | 15 | -0.070 | 0.805 |
| D7 | 14 | +0.380 | 0.180 |
| D30 | 15 | +0.400 | 0.140 |
| D90 | 15 | -0.265 | 0.341 |

**Placebo arm:**

| Timepoint | n | rho | p-value |
|---|---|---|---|
| D0 | 6 | +0.429 | 0.397 |
| D7 | 6 | +0.200 | 0.704 |
| D30 | 6 | +0.029 | 0.957 |
| D90 | 5 | -0.700 | 0.188 |

### Interpretation
No correlations reached statistical significance. The FMT arm shows a biologically interesting pattern â€” no correlation at baseline, positive trend at D7-D30 (both markers co-elevated in acute post-FMT destabilisation phase), then dissipation at D90. All results are hypothesis-generating given the underpowered sample sizes (feasibility trial).

---

## 8. Multi-Correlation Analysis

### Script written
File: `/scratch/prj/chmi_rbiome/project/scripts/MSc-Thesis/stats/PROFIT_correlations.R`

### What the script does
- Correlates 5 outcome variables against 13 faecal cytokines
- Runs separately for FMT and Placebo arms
- Runs separately for each timepoint (D0, D7, D30, D90)
- Applies Benjamini-Hochberg FDR correction across all tests
- Saves full results to CSV

### Outcome variables
E_faecalis_copies, E_coli_copies, MELD, Calprotectin, faecal_ammonia

### Cytokines tested
faecal_IL17A, IL17E, IL17F, IL21, IL22, IFNg, IL10, IL1b, IL6, TNFa, IL12, IL23, IL8

### Results
Results saved to: `/scratch/prj/chmi_rbiome/project/results/PROFIT_correlations_all.csv`

**Significant finding after FDR correction:**

| Arm | Timepoint | Outcome | Cytokine | n | rho | p_adjusted |
|---|---|---|---|---|---|---|
| Placebo | D0 | E_faecalis_copies | faecal_IL17F | 5 | 1.000 | 0.000 |

**Interpretation:** Perfect positive correlation between E. faecalis copy number and faecal IL-17F at baseline in the placebo group. Biologically plausible â€” E. faecalis is a pathobiont associated with Th17-driven inflammation. **Critical caveat: n=5. Requires validation in larger cohort.**

**All other correlations:** Non-significant after FDR correction. Expected given feasibility trial sample sizes (n=5-15 per arm per timepoint).

---

## 9. GitHub Commits This Session

| Commit message | What was saved |
|---|---|
| "Add Trimmomatic pipeline script for POIROT samples" | trimmomatic.sh |
| "Increase Trimmomatic time limit to 12hrs and add skip for completed samples" | Updated trimmomatic.sh |
| "Add PROFIT correlation analysis script" | PROFIT_correlations.R |
| "Add numeric conversion for correlation columns" | Updated PROFIT_correlations.R |

---

## 10. Outstanding Tasks and Blockers

| Task | Status | Notes |
|---|---|---|
| Trimmomatic â€” remaining 7 samples | â³ Running (job 34326136) | Check tomorrow morning |
| Bowtie2 script | âŒ Not started | Write after Trimmomatic completes |
| MetaPhlAn script | âŒ Not started | Write after Bowtie2 completes |
| ProxiMeta data transfer to HPC | âŒ Blocked | Awaiting Joseph/supervisor guidance |
| POIROT metadata | âŒ Blocked | No access yet |
| PROFIT Day7/Day30/Day90 CSVs | â³ Pending | Transferred to HPC, not yet analysed in R |
| Local R version check | â³ Pending | Not yet done |
| Supervisor confirmation on data sharing ethics | â³ Pending | |

---

## 11. Key Decisions Made This Session

1. **Plasma vs faecal ammonia** â€” switched from plasma (col 29) to faecal (col 78) for correlation with IL-22 as biologically more appropriate
2. **Spearman throughout** â€” both ammonia and IL-22 non-normal at all timepoints
3. **FDR correction applied** â€” Benjamini-Hochberg across all multi-correlations
4. **Novel focus confirmed** â€” correlations focused on E. coli, E. faecalis, MELD, Calprotectin vs cytokines â€” not reproducing pre-print cytokine vs cytokine work
5. **12-hour SLURM time limit** â€” increased from 4 hours after first job risked cancellation
6. **tmux session established** â€” named "thesis", prevents future disconnection losses

---

## 12. First Thing to Do Next Session

1. Check Trimmomatic job 34326136 â€” `squeue -u k25118483`
2. If complete, check log: `cat /scratch/prj/chmi_rbiome/project/logs/trimmomatic_34326136.log`
3. Verify all 10 samples output: `ls /scratch/prj/chmi_rbiome/project/results/trimmomatic/ | wc -l` (should show 40 files â€” 4 per sample Ã— 10 samples)
4. Write Bowtie2 script

---

*End of Session 2 â€” Next session: Verify Trimmomatic output, write Bowtie2 human read removal script*
