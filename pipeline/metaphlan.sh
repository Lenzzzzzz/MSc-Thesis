#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G

export PATH=/users/k25126777/.conda/envs/metaphlan/bin:$PATH

DB=/scratch/prj/chmi_rbiome/databases/Jan25_metaphlan_db
IN=/scratch/prj/chmi_rbiome/project/results/bowtie2
IN2=/scratch/users/k25118483/bowtie2
OUT=/scratch/users/k25118483

mkdir -p $OUT

for R1 in $IN/*_unmapped_R1.fastq.gz $IN2/*_unmapped_R1.fastq.gz; do
    SAMPLE=$(basename $R1 _unmapped_R1.fastq.gz)
    R2=$(dirname $R1)/${SAMPLE}_unmapped_R2.fastq.gz

    if [ -f "$OUT/${SAMPLE}_profile.tsv" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi

    echo "Processing: $SAMPLE"

    metaphlan $R1,$R2 \
        --input_type fastq \
        --db_dir $DB \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --nproc 1 \
        --output_file $OUT/${SAMPLE}_profile.tsv \
        --mapout $OUT/${SAMPLE}_bowtie2.bz2

    echo "Finished: $SAMPLE"
done

echo "All samples complete"
