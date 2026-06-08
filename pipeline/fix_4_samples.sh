#!/bin/bash
#SBATCH --job-name=fix_4_samples
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/fix_4_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/fix_4_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

module load trimmomatic/0.39-gcc-13.2.0
module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

RAW=/scratch/prj/chmi_rbiome/lindsey_data/kingco_007_0226
TRIM_OUT=/scratch/users/k25118483/trimmomatic
BT2_OUT=/scratch/users/k25118483/bowtie2
REF=/scratch/prj/chmi_rbiome/databases/GRCh38_noalt_as/GRCh38_noalt_as

mkdir -p $TRIM_OUT
mkdir -p $BT2_OUT

for SAMPLE in KINGCO-007-AKQ014B KINGCO-007-AKQ01B KINGCO-007-AKQ03B KINGCO-007-AKQ05A; do

    echo "=== Trimmomatic: $SAMPLE ==="
    trimmomatic PE -threads 4 \
        $RAW/${SAMPLE}_R1.fastq.gz \
        $RAW/${SAMPLE}_R2.fastq.gz \
        $TRIM_OUT/${SAMPLE}_R1_paired.fastq.gz \
        $TRIM_OUT/${SAMPLE}_R1_unpaired.fastq.gz \
        $TRIM_OUT/${SAMPLE}_R2_paired.fastq.gz \
        $TRIM_OUT/${SAMPLE}_R2_unpaired.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:50

    echo "=== Bowtie2: $SAMPLE ==="
    bowtie2 -x $REF \
        -1 $TRIM_OUT/${SAMPLE}_R1_paired.fastq.gz \
        -2 $TRIM_OUT/${SAMPLE}_R2_paired.fastq.gz \
        --threads 8 \
        --un-conc-gz $BT2_OUT/${SAMPLE}_unmapped_R%.fastq.gz \
        -S /dev/null

    echo "Finished: $SAMPLE"
done

echo "All 4 samples complete"
