#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

IN=/scratch/prj/chmi_rbiome/project/results/bowtie2
OUT=/scratch/prj/chmi_rbiome/project/results/metaphlan
DB=/scratch/prj/chmi_rbiome/databases/Jan25_metaphlan_db
INDEX=mpa_vJan25_CHOCOPhlAnSGB_202503

for R1 in $IN/*_unmapped_R1.fastq.gz; do
    SAMPLE=$(basename $R1 _unmapped_R1.fastq.gz)
    R2=$IN/${SAMPLE}_unmapped_R2.fastq.gz

 if [ -f "/scratch/users/k25118483/${SAMPLE}_profile.tsv" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi

    echo "Processing: $SAMPLE"

    metaphlan $R1,$R2 \
        --input_type fastq \
        --db_dir $DB \
        --index $INDEX \
        --nproc 8 \
        --mapout $OUT/${SAMPLE}.bowtie2.bz2 \
        -o $OUT/${SAMPLE}_profile.tsv

cp $OUT/${SAMPLE}_profile.tsv /scratch/users/k25118483/${SAMPLE}_profile.tsv
    echo "Finished: $SAMPLE"
done

echo "All samples complete"
