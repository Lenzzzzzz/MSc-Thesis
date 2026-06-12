# Computational Lab Notebook
## MSc Thesis — Investigating the 'Immune-AMR' Axis
**Author:** Eleni Andria Kalopedis
**Supervisor:** Dr. Lindsey Ann Edwards, King's College London
**Session:** 6
**Date:** 10–12 June 2026
**Platform:** CREATE HPC (King's College London) + Chromebook/Linux local

---

## Session Overview

This session completed the POIROT MetaPhlAn analysis, performed a comprehensive ProxiMeta ARG analysis using the PME bulk resistance export, conducted bin-level metabolic investigations for three paired patients (AKQ003, AKQ005, AKQ014), ran four statistical tests, and generated 12 new figures. The most significant findings include the identification of an ExPEC virulence profile in AKQ003B, characterisation of E. coli metabolic fitness and HGT capacity, and three distinct antibiotic response phenotypes across the cohort.

---

## 1. ProxiMeta Data — Files Downloaded and Transferred

### Standard ProxiMeta Platform
Individual tar.gz archives downloaded per sample:
- `amr_files.tar.gz` — contains full AMR data including matrix files
- `plasmid_files.tar.gz` — plasmid host associations
- `viral_files.tar.gz` — viral MAG data
- `microbial_genome_files.tar.gz` — proximiteta_report.tsv (MAG quality)

Metabolism files downloaded for AKQ003A, AKQ003B, AKQ005A, AKQ005B, AKQ014A, AKQ014B and transferred to:
```
/scratch/prj/chmi_rbiome/project/metadata/PROXIMETA_POIROT/
AKQ_003_A_Metabolism_Files/metabolism_results/cluster_modules_matrix_complete.tsv
AKQ_003_B_Metabolism_Files/metabolism_results/cluster_modules_matrix_complete.tsv
AKQ_005_A_Metabolism_Files/metabolism_results/cluster_modules_matrix_complete.tsv
AKQ_005B_Metabolism_Files/metabolism_results/cluster_modules_matrix_complete.tsv
AKQ_014_A_Metabolism_Files/metabolism_results/cluster_modules_matrix_complete.tsv
AKQ_014_B_Metabolism_Files/metabolism_results/cluster_modules_matrix_complete.tsv
```

### ProxiMeta Explorer (PME)
Key files transferred to HPC:
- `POIROTallsamples.csv` — species abundance + MAG quality, 9432 rows, all POIROT samples
- `POIROTallsamplesresistance.csv` — PME resistance export (amr_composition equivalent), 557 rows, 7 samples
- Individual sample CSVs: AKQ002A, AKQ003A/B, AKQ005A/B, AKQ014A/B.csv
- `PROFIT.csv` — sample catalogue only (not usable for analysis alone)

**Key clarification:** POIROTallsamplesresistance.csv is the primary ARG import file per M1T2 framework. It replaces the three-file join and includes full species-level taxonomy, Contig TPM, Genome TPM, Host Completeness/Contamination, and Gene contig type (plasmid/genomic/viral) pre-computed.

---

## 2. PME Resistance Analysis — Key Findings

### File structure confirmed
- 557 rows, 23 columns
- 7 samples: AKQ_014_A, AKQ_002_A, AKQ_014_B, AKQ_003_B, AKQ_003_A, AKQ_005B, AKQ_005_A
- Element types: AMR (304), STRESS (129), VIRULENCE (124)
- Gene contig types: genomic (271), plasmid (101), viral (5)

### Viral ARG revision
**Correction to earlier analysis:** 1 genuine viral AMR gene detected — erm(B) in AKQ_014_A on a viral contig in Hominilimicola fabiformis (91.95% completeness). The other 4 viral entries are VIRULENCE genes, not AMR. The earlier finding of "zero viral ARGs" from individual files was therefore incorrect. Plasmids remain the dominant MGE pathway but are not the exclusive one.

### Top ARG-carrying species (PME species-resolved)
| Rank | Species | Total ARGs | n_samples |
|---|---|---|---|
| 1 | NoHost | 100 | 7 |
| 2 | Escherichia coli | 15 | 5 |
| 3 | Alistipes putredinis | 10 | 5 |
| 4 | Bacteroides uniformis | 9 | 5 |

---

## 3. Bin-Level Deep-Dive Analysis

### AKQ005A — E. coli bin_10
- Species: Escherichia coli
- Completeness: 73.63%, Contamination: 3.71% (MQ)
- Genome TPM: 2696.6
- AMR genes: blaEC (chromosomal beta-lactamase), emrD (efflux)
- STRESS genes: Complete arsenic operon (arsABCDR), acid resistance (ariR plasmid, asr genomic), iron efflux (fieF), copper tolerance
- VIRULENCE genes: cdtB (cytolethal distending toxin), espX1 (PLASMID — type III secretion effector), lpfA/lpfA-O113 (colonisation fimbriae)
- Metabolic pathways complete: ~115 including glycolysis, PPP, TCA, amino acid biosynthesis
- Type IV secretion: COMPLETE — conjugative plasmid transfer capable
- Multiple complete efflux systems: MdtABC, MdlAB, EmrAB, MdtIJ
- EvgS-EvgA (acid tolerance): COMPLETE

### AKQ005 T1 vs T2 comparison
| Feature | T1 (AKQ005A) | T2 (AKQ005B) |
|---|---|---|
| AMR genes | 44 | 50 |
| STRESS genes | 12 | 19 |
| VIRULENCE genes | 29 | 16 |
| Plasmid ARG fraction | 25% | 12% |
| E. coli present | Yes (TPM 2696) | **No — eliminated** |
| Dominant organism | Alistipes putredinis | Alistipes putredinis |
| Dominant org. complete pathways | 115 (E. coli) | 39 (Alistipes) |
| Type IV secretion | Complete | Absent |

**AKQ005 phenotype:** Pathobiont elimination. Antibiotic eliminated HGT-competent E. coli but caused massive community collapse (174 extinctions). Post-antibiotic community dominated by simpler, less resistance-capable commensal.

### AKQ014 T1 vs T2 comparison
| Feature | T1 (AKQ014A) | T2 (AKQ014B) |
|---|---|---|
| AMR genes | 51 | 52 |
| STRESS genes | 24 | **58 (2.4×↑)** |
| VIRULENCE genes | 22 | 21 |
| E. coli TPM | 2456.9 | 1354.8 (impaired) |
| E. coli completeness | 66.2% | 4.17% (fragment) |
| E. coli metabolic pathways | 93 | 17 |
| ESKAPE colonisation | None | K. pneumoniae, E. faecalis |

**AKQ014 phenotype:** E. coli metabolic collapse + simultaneous ESKAPE invasion. The tripling of STRESS genes indicates community-wide antibiotic stress response.

**Note on species discrepancy:** MetaPhlAn detected Enterococcus faecium; PME resolves to Enterococcus faecalis. Methodological difference documented.

### AKQ003 T1 vs T2 comparison
| Feature | T1 (AKQ003A) | T2 (AKQ003B) |
|---|---|---|
| AMR genes | 32 | 31 |
| STRESS genes | 6 | 5 |
| VIRULENCE genes | 2 | **14 (7×↑)** |
| E. coli present | No | **Yes — new colonisation (TPM 3286)** |
| Dominant T1 organism | Ruminococcus bromii_B | — |
| Dominant T2 organism | — | Escherichia coli (2 bins) |

**AKQ003B E. coli virulence profile — ExPEC:**
- ibeA — brain endothelium invasion (meningitis/sepsis)
- sfaF/sfaS — S fimbriae (UTI/meningitis pathogenicity)
- iroB/C/D/E/N — salmochelin siderophore (iron acquisition in blood)
- iss — increased serum survival (complement resistance)
- ybtP/Q — yersiniabactin (iron acquisition)
- vactox — vacuolating cytotoxin
- fdeC, sslE — mucosal adhesion and colonisation

**This is an ExPEC (Extraintestinal Pathogenic E. coli) virulence profile** associated with meningitis, UTIs and sepsis. It colonised post-antibiotic in a surgical patient. This is the most clinically alarming finding in the POIROT cohort.

AKQ003B E. coli metabolic capacity: 57 complete pathways, efflux systems (AcrAB-TolC, EmrAB) complete, no Type IV secretion.

**AKQ003 phenotype:** De novo ExPEC colonisation post-antibiotic with 7-fold virulence expansion.

---

## 4. Three Antibiotic Response Phenotypes

| Phenotype | Patient | Key feature |
|---|---|---|
| Pathobiont elimination | AKQ005 | E. coli eliminated, 174 extinctions, HGT capacity lost |
| E. coli collapse + ESKAPE invasion | AKQ014 | E. coli impaired, K. pneumoniae + E. faecalis colonise, STRESS 2.4×↑ |
| De novo ExPEC colonisation | AKQ003 | ExPEC E. coli new dominant organism, virulence 7×↑ |

---

## 5. Statistical Tests Run

### 1. Spearman: Genome TPM vs n_ARGs
- rho = 0.123, p = 0.243, n = 92 bins
- **Interpretation:** No significant correlation. ARG carriage is independent of organism abundance. ARG surveillance cannot rely on abundance profiling alone.

### 2. Fisher's Exact Test: ESKAPE colonisation T1 vs T2
- T1: 0/5 patients; T2: 3/5 patients (incomplete — 2 pending)
- p = 0.1667 (not yet significant)
- **Critical note:** When AKQ001B and AKQ002B are available, if neither carries ESKAPE pathogens, p = 0.0476 (significant). RE-RUN when all 5 T2 samples available.

### 3. Descriptive statistics: ARG burden
- T1 mean: 42.8 (95% CI: 30.2–55.3), n=4 samples
- T2 mean: 44.3 (95% CI: 15.5–73.1), n=3 samples
- MGE fraction T1: 17.0%; T2: 15.6%
- Wide T2 CI reflects incomplete n

### 4. Chi-square: plasmid vs genomic ARG distribution
- chi = 0, p = 1.0
- T1 plasmid fraction: 24.8%; T2: 25.6%
- **Note:** Pooled analysis loses within-patient signal. Correct test is paired Wilcoxon on MGE fractions — pending n=5 complete pairs.

### Results saved to:
`/scratch/prj/chmi_rbiome/project/results/POIROT_statistical_tests.csv`

---

## 6. Figures Generated This Session

All saved to `/scratch/prj/chmi_rbiome/project/figures/`

| Figure | Description |
|---|---|
| AKQ005_element_type_comparison.pdf | AMR/STRESS/VIRULENCE T1 vs T2 |
| AKQ005_species_ARG_comparison.pdf | Top ARG-carrying species T1 vs T2 |
| AKQ005_MGE_fraction.pdf | Plasmid vs chromosomal ARG fraction |
| AKQ005_metabolic_comparison.pdf | Dominant organism pathway completeness |
| AKQ003_element_type.pdf | AMR/STRESS/VIRULENCE T1 vs T2 |
| AKQ003_species_ARG.pdf | Top ARG-carrying species — E. coli colonisation |
| AKQ003_MGE_fraction.pdf | Plasmid vs chromosomal ARG fraction |
| AKQ003_metabolic.pdf | Ruminococcus vs E. coli pathway completeness |
| AKQ014_element_type.pdf | AMR/STRESS/VIRULENCE — STRESS 2.4×↑ |
| AKQ014_species_ARG.pdf | Top ARG-carrying species T1 vs T2 |
| AKQ014_MGE_fraction.pdf | Plasmid vs chromosomal ARG fraction |
| AKQ014_metabolic.pdf | E. coli metabolic collapse T1→T2 |

---

## 7. Analysis Status — Outstanding

### Immediate (pending AKQ001B and AKQ002B):
- [ ] Re-run Fisher's exact test (ESKAPE) — will likely reach p<0.05
- [ ] Wilcoxon signed-rank: ARG burden T1 vs T2 (n=5 pairs)
- [ ] Wilcoxon signed-rank: MGE fraction T1 vs T2
- [ ] Wilcoxon signed-rank: VIRULENCE count T1 vs T2
- [ ] AKQ001 and AKQ002 bin-level analysis

### Downloads still needed (standard ProxiMeta platform):
- [ ] plasmid_files.tar.gz per sample (plasmid host associations)
- [ ] viral_files.tar.gz per sample (formal viral ARG validation)
- [ ] proximiteta_report.tsv per sample (MAG quality — supervisor requested)

### Analysis pending:
- [ ] Spearman: species TPM vs ARG burden (ecological overlay)
- [ ] Integration with POIROT immune/inflammatory data
- [ ] Cross-patient summary figure (all 5 patients)
- [ ] PROFIT PME bulk export and resistome analysis
- [ ] PROFIT cluster × ARG correlation

### Writing:
- [ ] Methods section — computational pipeline
- [ ] Results 3.2.4 update — ExPEC finding (AKQ003B)
- [ ] Results 3.2.5 update — bin-level metabolic narrative
- [ ] Discussion — three antibiotic phenotypes
- [ ] Discussion — ExPEC colonisation clinical implications

---

## 8. First Things To Do Next Session

```bash
# Check if AKQ001B and AKQ002B have finished uploading
# Log into ProxiMeta platform and check Organisation Samples

# If available, download AMR files for AKQ001B and AKQ002B
# Then re-run Fisher's exact test

# Update results section with ExPEC and metabolic findings
```

---

*End of Session 6 — Next session: Complete paired statistical analyses once AKQ001B and AKQ002B available, update results section with bin-level findings, begin PROFIT resistome analysis*
