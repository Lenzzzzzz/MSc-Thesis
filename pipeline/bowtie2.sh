#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/bowtie2_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/bowtie2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

TRIMMED=/scratch/prj/chmi_rbiome/project/results/trimmomatic
OUT=/scratch/users/k25118483
DB=/scratch/prj/chmi_rbiome/databases/GRCh38_noalt_as/GRCh38_noalt_as

for SAMPLE in KINGCO-007-AKQ014B KINGCO-007-AKQ01B KINGCO-007-AKQ03B KINGCO-007-AKQ05A; do
    R1=$TRIMMED/${SAMPLE}_R1_paired.fastq.gz
    R2=$TRIMMED/${SAMPLE}_R2_paired.fastq.gz

    if [ -f "$OUT/${SAMPLE}_unmapped_R1.fastq.gz" ]; then
        echo "Skipping $SAMPLE - already done"
        continue
    fi

    echo "Processing: $SAMPLE"

    bowtie2 -x $DB \
        -1 $R1 \
        -2 $R2 \
        --threads 8 \
        --sensitive \
        --un-conc-gz $OUT/${SAMPLE}_unmapped_R%.fastq.gz \
        -S /dev/null \
        2>> /scratch/prj/chmi_rbiome/project/logs/bowtie2_stats.log

    echo "Finished: $SAMPLE"
done

echo "All samples complete"
EOF
