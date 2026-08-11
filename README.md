# The Immune-AMR Axis: How Antibiotic Exposure and Chronic Inflammation Shape Gut Resistome Dynamics

**Author:** Eleni Andria Kalopedis  
**Supervisor:** Dr Lindsey Ann Edwards, King's College London  
**Programme:** MSc Microbiome in Health and Disease, 2026  

---

## Overview

This repository contains all analysis scripts for the MSc thesis investigating
how the host inflammatory microenvironment shapes gut resistome composition,
mobile ARG dynamics, and immune-resistome associations across two clinical cohorts:

- **POIROT** — patients with orthopaedic trauma who developed drug-resistant
  wound infections requiring antibiotic treatment (acute perturbation model, n=5)
- **PROFIT** — randomised FMT trial in patients with liver cirrhosis
  (chronic inflammatory model, n=21; Woodhouse et al., 2019)

Hi-C proximity ligation (ProxiMeta, Phase Genomics) was used to attribute ARGs
to specific bacterial hosts and resolve chromosomal versus plasmid encoding.

---

## Repository Structure

**pipeline/** — HPC bash scripts for quality trimming (Trimmomatic),
host depletion (Bowtie2), and taxonomic profiling (MetaPhlAn 4)

**stats/** — R analysis scripts:
- `POIROT_ARG_analysis.R` — Resistome characterisation
- `POIROT_metaphlan_analysis.R` — Diversity metrics and PERMANOVA
- `PROFIT_ARG_analysis.R` — Baseline resistome characterisation
- `PROFIT_longitudinal_analysis.R` — MGE fraction and FMT vs placebo
- `PROFIT_correlations.R` — Spearman immune-ARG correlations
- `PROFIT_immune_correlations.R` — MSD ECL and LE assay integration
- `PROFIT_clinical_analysis.R` — GGT, CRP, ammonia analysis
- `PROFIT_patient_metabolic_analysis.R` — Patient-level metabolic capacity

**figures/** — Figure generation scripts (ggplot2)

**lab_notebook/** — Computational lab notebook documenting analysis sessions,
pipeline decisions, and troubleshooting steps

**thecleanup/** — Metadata preprocessing scripts; archived after integration
into main analysis pipeline

---

## Data Availability

Patient data are not included in this repository in accordance with
Material Transfer Agreement obligations and KCL data governance policy.
POIROT samples were provided by Professor Ana Valdes and Dr Afroditi Kouraki,
Nottingham University Hospitals NHS Trust, under a Material Transfer Agreement.
PROFIT trial data are described in Woodhouse et al. (2019) BMJ Open.

All analysis was performed on the KCL CREATE high-performance computing cluster.

---

## Dependencies

R v4.5.2. Full package list in Supplementary Table 1C of the thesis.  
Key packages: vegan, ggplot2, dplyr, tidyr, readr, pheatmap, ggrepel, patchwork.
