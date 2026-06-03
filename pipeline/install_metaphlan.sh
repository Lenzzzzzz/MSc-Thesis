#!/bin/bash
#SBATCH --job-name=install_mpa
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/install_mpa_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/install_mpa_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

module load anaconda3/2022.10-gcc-13.2.0
source $(conda info --base)/etc/profile.d/conda.sh

conda env remove -n mpa -y 2>/dev/null

conda create --name mpa -c conda-forge -c bioconda python=3.7 metaphlan -y

conda activate mpa
metaphlan --version
echo "Installation complete"
which metaphlan
EOF
