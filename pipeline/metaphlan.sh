#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G

module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

DB=/scratch/prj/chmi_rbiome/databases/Jan25_metaphlan_db
OUT=/scratch/users/k25118483

mkdir -p $OUT

# Process 6 samples from Joseph's bowtie2 output
for SAMPLE in KINGCO-007-AKQ014A KINGCO-007-AKQ01A KINGCO-007-AKQ02A KINGCO-007-AKQ02B KINGCO-007-AKQ03A KINGCO-007-AKQ05B; do
    R1=/scratch/prj/chmi_rbiome/project/results/bowtie2/${SAMPLE}_unmapped_R1.fastq.gz
    R2=/scratch/prj/chmi_rbiome/project/results/bowtie2/${SAMPLE}_unmapped_R2.fastq.gz

    if [ -f "$OUT/${SAMPLE}_profile.tsv" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi

    echo "Processing: $SAMPLE"
    metaphlan $R1,$R2 \
        --input_type fastq \
        --db_dir $DB \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --nproc 4 \
        --output_file $OUT/${SAMPLE}_profile.tsv \
        --mapout $OUT/${SAMPLE}_bowtie2.bz2
    echo "Finished: $SAMPLE"
done

# Process 4 samples from our fixed bowtie2 output
for SAMPLE in KINGCO-007-AKQ014B KINGCO-007-AKQ01B KINGCO-007-AKQ03B KINGCO-007-AKQ05A; do
    R1=/scratch/users/k25118483/bowtie2/${SAMPLE}_unmapped_R1.fastq.gz
    R2=/scratch/users/k25118483/bowtie2/${SAMPLE}_unmapped_R2.fastq.gz

    if [ -f "$OUT/${SAMPLE}_profile.tsv" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi

    echo "Processing: $SAMPLE"
    metaphlan $R1,$R2 \
        --input_type fastq \
        --db_dir $DB \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --nproc 4 \
        --output_file $OUT/${SAMPLE}_profile.tsv \
        --mapout $OUT/${SAMPLE}_bowtie2.bz2
    echo "Finished: $SAMPLE"
done

echo "All samples complete"
