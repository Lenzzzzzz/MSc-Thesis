# Computational Lab Notebook
## MSc Thesis — Investigating the 'Immune-AMR' Axis
**Author:** Eleni Andria Kalopedis
**Supervisor:** Dr. Lindsey Ann Edwards, King's College London
**Session:** 4
**Date:** 2 June 2026
**Platform:** CREATE HPC (King's College London) + Chromebook/Linux local

---

## Session Overview

This session focused on completing the Trimmomatic pipeline, resolving the Bowtie2 memory failure, installing MetaPhlAn correctly, and troubleshooting the MetaPhlAn database path. Key findings included a data labelling inconsistency (AKQ014B/AKQ041B), MD5 checksum verification of all raw files, and identification that Joseph's Bowtie2 output files are more complete and should be used for MetaPhlAn input.

---

## 1. Trimmomatic — Final Completion

### Background
Previous session left AKQ03B and AKQ05A/AKQ05B incomplete due to a node-specific filesystem issue (erc-hpc-comp016).

### Resolution
Resubmitted to a different node (erc-hpc-comp003). AKQ03B completed successfully — confirming the issue was node-specific, not data corruption.

### MD5 Checksum Verification
All 20 raw sequencing files verified against their MD5 checksums:

```bash
for sample in AKQ014A AKQ014B AKQ01A AKQ01B AKQ02A AKQ02B AKQ03A AKQ03B AKQ05A AKQ05B; do
    md5sum -c md5/KINGCO-007-${sample}_R1.md5 2>/dev/null
    md5sum -c md5/KINGCO-007-${sample}_R2.md5 2>/dev/null
done
```

**Results:** All files returned OK except AKQ041B — which FAILED because the MD5 file references AKQ041B but the actual file is named AKQ014B.

### AKQ014B / AKQ041B Labelling Inconsistency

**Finding:** The MD5 file for AKQ014B internally references the filename `KINGCO-007-AKQ041B` (digits 14 transposed to 41). The same transposition appears in the pre-existing merged MetaPhlAn abundance table column headers.

**Impact:** When merging pipeline outputs with the existing MetaPhlAn table, AKQ014B and AKQ041B must be treated as the same sample. A manual rename or mapping will be required at the merge step.

**Action:** Flagged to Joseph via WhatsApp. Awaiting confirmation.

**Documented in:** `lab_notebook/KNOWN_ISSUES.md`

### Final Trimmomatic Statistics — All 10 Samples

| Sample | Input Pairs | Both Surviving | % Retained | Dropped |
|---|---|---|---|---|
| AKQ014A | 87,108,174 | 85,516,222 | 98.17% | 0.46% |
| AKQ014B | 78,029,094 | 76,620,986 | 98.20% | 0.50% |
| AKQ01A | 28,707,007 | 27,908,711 | 97.22% | 0.50% |
| AKQ01B | 88,894,007 | 87,057,017 | 97.93% | 0.83% |
| AKQ02A | 88,058,796 | 86,273,588 | 97.97% | 0.62% |
| AKQ02B | 114,872,412 | 111,285,577 | 96.88% | 1.75% |
| AKQ03A | 80,106,982 | 78,623,517 | 98.15% | 0.63% |
| AKQ03B | 66,420,740 | 65,187,189 | 98.14% | 0.53% |
| AKQ05A | 108,317,590 | 105,833,497 | 97.71% | 1.00% |
| AKQ05B | 74,319,881 | 72,901,898 | 98.09% | 0.53% |

**Mean retention: 97.95% | Range: 96.88%–98.20%**

**Thesis methods statement:**
*"Raw sequencing reads were quality trimmed using Trimmomatic v0.39 in paired-end mode (ILLUMINACLIP: TruSeq3-PE adapters 2:30:10; LEADING:3; SLIDINGWINDOW:4:15; TRAILING:3; MINLEN:36). A mean of 97.95% of read pairs were retained across all 10 samples (range: 96.88–98.20%)."*

---

## 2. Bowtie2 — Memory Failure and Resolution

### What happened
The Bowtie2 job (34386849) failed with exit code 137 — killed by the system due to insufficient memory. The `--mem=32G` allocation was insufficient for aligning large metagenomic samples (up to 114 million read pairs) against the full GRCh38 human genome.

### Error
```
(ERR): bowtie2-align exited with value 137
Killed
```

### Root cause
Memory parameter set to 32G — too low for this data size. Correct allocation is 64G minimum.

### Note on responsibility
The 32G memory allocation was set incorrectly in the initial script. This is documented as a calibration error. The script has been updated to 64G for future runs.

### Resolution
Joseph Falconer (k25126777) independently ran Bowtie2 on the same samples with a more comprehensive approach. His output files (`_unmapped_R1/R2.fastq.gz`) are present in the bowtie2 results folder and were used for MetaPhlAn input.

### File size comparison
| File type | AKQ014A R1 size | Notes |
|---|---|---|
| Our `_microbial1_1.fastq.gz` | 271 MB | Only concordantly unaligned pairs |
| Joseph's `_unmapped_R1.fastq.gz` | 858 MB | More comprehensive host removal |

### Decision
Joseph's `_unmapped` files used for MetaPhlAn input. This is appropriate collaborative practice and does not affect thesis validity.

### Thesis methods statement
*"Following quality trimming, human host reads were removed by alignment to the GRCh38 reference genome using Bowtie2 v2.5.1. Unaligned reads representing the microbial fraction were retained for downstream taxonomic profiling."*

---

## 3. MetaPhlAn Installation

### Installation attempts

**Attempt 1:** `conda install -c bioconda metaphlan` on login node — killed by HPC (insufficient memory for conda solve on login node)

**Attempt 2:** SLURM installation job using `--prefix /scratch/prj/chmi_rbiome/project/mpa_env` — failed silently. Conda environment not created at specified path. pip-installed MetaPhlAn v4.2.4 at `/users/k25118483/.local/bin/metaphlan` was incorrectly reported as successful.

**Attempt 3:** SLURM installation job using exact wiki command:
```bash
conda create --name mpa -c conda-forge -c bioconda python=3.7 metaphlan -y
```
Job ID: 34486868 — running at end of session.

### Key lesson
The pip-installed MetaPhlAn (v4.2.4) has different flag names from the conda-installed version. Specifically:
- pip version does NOT support `--bowtie2db` or `--bowtie2out`
- conda version DOES support these flags as documented in the wiki

### MetaPhlAn flag history and corrections

| Attempt | Flags used | Result |
|---|---|---|
| 1 | `--bowtie2db`, `--bowtie2out`, `--unclassified_estimation` | Failed — pip version doesn't support these |
| 2 | `--database`, `--mapout` | Failed — wrong flag names |
| 3 | `DEFAULT_DB_FOLDER` env variable | Failed — MetaPhlAn ignored it, tried to download database to home directory |
| 4 | `--offline`, `--index $DB/$INDEX` | Failed — MetaPhlAn still looked in default directory |
| 5 | Conda env + `--bowtie2db`, `--bowtie2out` | Pending — correct approach per wiki and Joseph |

### Disk quota issue
MetaPhlAn repeatedly tried to download the database (30GB+) into `/users/k25118483/` which only has 50GB quota. This filled the home directory. Resolved by:
1. Deleting partial database downloads
2. Using `--offline` flag (later)
3. Using conda environment which installs in scratch

### Joseph's guidance (WhatsApp, 1 June 2026)
*"You should use a conda environment (not your base env) to install and run metaphlan. If you need to install metaphlan use the conda install function or pip install like on the github that Lindsey shared above. Then you don't need to download the database again, just use the -x tag to tell metaphlan where to find the database."*

### Correct MetaPhlAn command (per wiki + Joseph)
```bash
metaphlan $R1,$R2 \
    --bowtie2db /scratch/prj/chmi_rbiome/databases/Jan25_metaphlan_db \
    --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
    --input_type fastq \
    --nproc 8 \
    --bowtie2out $OUT/${SAMPLE}.bowtie2.bz2 \
    -o $OUT/${SAMPLE}_profile.tsv
```

---

## 4. FastQC — Standalone Script

FastQC standalone script not yet written. To be completed next session after MetaPhlAn installation confirmed working.

---

## 5. GitHub Commits This Session

| Commit message | What was saved |
|---|---|
| "Fix bowtie2 output filename bug in un-conc-gz parameter" | bowtie2.sh |
| "Update Trimmomatic stats with final 3 samples and node issue documentation" | trimmomatic_stats_formatted.txt |
| "Add MetaPhlAn taxonomic profiling script for POIROT samples" | metaphlan.sh |
| "Add Bowtie2 module load to MetaPhlAn script" | metaphlan.sh |
| "Fix MetaPhlAn flags for version 4.2.4" | metaphlan.sh |
| "Fix MetaPhlAn database path using DEFAULT_DB_FOLDER variable" | metaphlan.sh |
| "Fix MetaPhlAn to use full database path and offline mode" | metaphlan.sh |
| "Fix MetaPhlAn to use conda environment and correct bowtie2db flag per wiki" | metaphlan.sh |
| "Fix MetaPhlAn conda installation script using exact wiki command" | install_metaphlan.sh |

---

## 6. Known Issues (updated)

| Issue | Status |
|---|---|
| AKQ014B/AKQ041B labelling inconsistency | Flagged to Joseph — awaiting confirmation |
| AKQ03B node-specific processing failure (erc-hpc-comp016) | Resolved — completed on erc-hpc-comp003 |
| Bowtie2 memory failure (32G insufficient) | Documented — Joseph's files used instead |
| MetaPhlAn conda installation | Pending — job 34486868 running |

---

## 7. Outstanding Tasks

| Task | Status | Notes |
|---|---|---|
| MetaPhlAn installation | ⏳ Running (job 34486868) | Check log for completion |
| MetaPhlAn profiling — all 10 samples | ❌ Pending | Submit after installation confirmed |
| FastQC standalone script | ❌ Not written | Write next session |
| Merge MetaPhlAn tables | ❌ Pending | After all profiles generated |
| POIROT metadata | ❌ Blocked | No access yet |
| ProxiMeta data transfer | ❌ Blocked | Awaiting supervisor/Joseph guidance |

---

## 8. First Things To Do Next Session

1. Check MetaPhlAn installation: `cat /scratch/prj/chmi_rbiome/project/logs/install_mpa_34486868.log | tail -10`
2. If complete, find MetaPhlAn binary: `which metaphlan` after `conda activate mpa`
3. Update `metaphlan.sh` to use full conda binary path
4. Submit MetaPhlAn job
5. Write standalone FastQC script

---

*End of Session 4 — Next session: Confirm MetaPhlAn installation, run MetaPhlAn on all 10 POIROT samples, write FastQC script*
