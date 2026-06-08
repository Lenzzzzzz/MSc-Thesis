#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/metaphlan_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G

# Initialize conda and activate YOUR original working environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metaphlan_env

DB=/scratch/prj/chmi_rbiome/databases/Jan25_metaphlan_db
IN=/scratch/prj/chmi_rbiome/project/results/bowtie2
IN2=/scratch/users/k25118483/bowtie2
OUT=/scratch/users/k25118483

mkdir -p $OUT

for R1 in $IN/*_unmapped_R1.fastq.gz $IN2/*_unmapped_R1.fastq.gz; do
    # Guard clause to skip if no files match a directory yet
    [ -e "$R1" ] || continue

    SAMPLE=$(basename $R1 _unmapped_R1.fastq.gz)
    R2=$(dirname $R1)/${SAMPLE}_unmapped_R2.fastq.gz
    
    # Skip condition: if your perfect profile already exists, do not rerun it!
    if [ -f "$OUT/${SAMPLE}_profile.tsv" ]; then
        echo "Skipping $SAMPLE - already perfectly processed"
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
