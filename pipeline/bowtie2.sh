#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/bowtie2_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/bowtie2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Load modules
module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

# Define paths
TRIMMED=/scratch/prj/chmi_rbiome/project/results/trimmomatic
OUT=/scratch/prj/chmi_rbiome/project/results/bowtie2
DB=/scratch/prj/chmi_rbiome/databases/GRCh38_noalt_as/GRCh38_noalt_as

# Loop through all samples
for R1 in $TRIMMED/*_R1_paired.fastq.gz; do
    SAMPLE=$(basename $R1 _R1_paired.fastq.gz)
    R2=$TRIMMED/${SAMPLE}_R2_paired.fastq.gz

    # Skip if already completed
    if [ -f "$OUT/${SAMPLE}_microbial_1.fastq.gz" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi

    echo "Processing: $SAMPLE"

    # Align to human genome and keep unaligned (microbial) reads
    bowtie2 -x $DB \
        -1 $R1 \
        -2 $R2 \
        --threads 8 \
        --sensitive \
        --un-conc-gz $OUT/${SAMPLE}_microbial1_%.fastq.gz \
        -S /dev/null \
        2>> /scratch/prj/chmi_rbiome/project/logs/bowtie2_stats.log

    echo "Finished: $SAMPLE"
done

echo "All samples complete - host reads removed"
