  /scratch/prj/chmi_rbiome/project/scripts/MSc-Thesis/pipeline/install_metaphlan.sh          
#!/bin/bash
#SBATCH --job-name=install_mpa
#SBATCH --output=/scratch/prj/chmi_rbiome/project/logs/install_mpa_%j.log
#SBATCH --error=/scratch/prj/chmi_rbiome/project/logs/install_mpa_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Load anaconda
module load anaconda3/2022.10-gcc-13.2.0

# Initialise conda for use in bash script
source $(conda info --base)/etc/profile.d/conda.sh

# Create MetaPhlAn environment in scratch space
conda create -y \
    --prefix /scratch/prj/chmi_rbiome/project/mpa_env \
    -c conda-forge \
    -c bioconda \
    python=3.9 \
    metaphlan

echo "MetaPhlAn installation complete"

# Test it works
conda activate /scratch/prj/chmi_rbiome/project/mpa_env
metaphlan --version
