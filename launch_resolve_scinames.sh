#!/bin/sh
#SBATCH -c 1
#SBATCH --mem 100G
#SBATCH -t 1-00:00:00
#SBATCH --job-name resolve
#SBATCH -p ulow

directory=$1
datafile=$2

echo "#### sourcing bashrc"
source ~/.bashrc
echo "#### activating conda env"
conda activate r-taxonomy-tools

cd $directory

echo "#### launching R calculation ########"
Rscript /home/lescailab/CODE/hgtseq_tools/resolve_scientific_names.R $datafile