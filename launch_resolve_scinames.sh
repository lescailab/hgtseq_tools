#!/bin/sh
#SBATCH -c 1
#SBATCH --mem 100G
#SBATCH -t 1-00:00:00
#SBATCH --job-name resolve
#SBATCH -p ulow

echo "#### sourcing bashrc"
source ~/.bashrc
echo "#### activating conda env"
conda activate r-taxonomy-tools

cd /home/lescailab/COLLABS/malacrida/glossina_hgtseq_run02

echo "#### launching R calculation ########"
Rscript resolve_scientific_names.R