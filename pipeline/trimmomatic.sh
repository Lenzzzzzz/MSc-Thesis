#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/trimmomatic_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/trimmomatic_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# Load Trimmomatic module
module load trimmomatic/0.39-gcc-13.2.0

# Define paths
RAW=/scratch/prj/chmi_rbiome/lindsey_data/kingco_007_0226
OUT=/scratch/prj/chmi_rbiome/project/results/trimmomatic
ADAPTERS=/scratch/prj/chmi_rbiome/databases/trimmomatic/TruSeq3-PE.fa

# Loop through all samples
for R1 in $RAW/*_R1.fastq.gz; do
    SAMPLE=$(basename $R1 _R1.fastq.gz)
    R2=$RAW/${SAMPLE}_R2.fastq.gz

if [ -f "$OUT/${SAMPLE}_R1_paired.fastq.gz" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi
    trimmomatic PE -threads 4 -phred33 \
        $R1 $R2 \
        $OUT/${SAMPLE}_R1_paired.fastq.gz \
        $OUT/${SAMPLE}_R1_unpaired.fastq.gz \
        $OUT/${SAMPLE}_R2_paired.fastq.gz \
        $OUT/${SAMPLE}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 \
        SLIDINGWINDOW:4:15 \
        TRAILING:3 \
        MINLEN:36
 2>&1 | tee -a /scratch/prj/chmi_rbiome/project/logs/trimmomatic_stats.log
    echo "Finished: $SAMPLE"
# ============================================================
# RUN FASTQC ON TRIMMED OUTPUT
# ============================================================

module load fastqc/0.12.1-gcc-13.2.0

# Create FastQC output directory
mkdir -p /scratch/prj/chmi_rbiome/project/results/fastqc

echo "Running FastQC on trimmed paired reads..."

fastqc \
    $OUT/*_R1_paired.fastq.gz \
    $OUT/*_R2_paired.fastq.gz \
    --outdir /scratch/prj/chmi_rbiome/project/results/fastqc \
    --threads 4

echo "FastQC complete"
done

echo "All samples complete"

