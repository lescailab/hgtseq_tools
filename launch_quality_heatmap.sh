#!/bin/sh
#SBATCH -c 16
#SBATCH --mem 100G
#SBATCH -t 2-00:00:00
#SBATCH --job-name hgtseq_quality
#SBATCH -p ulow

directory=$1
datafile=$2
parallel=$3

echo "#### sourcing bashrc"
source ~/.bashrc
echo "#### activating conda env"
conda activate r-taxonomy-tools

cd $directory

echo "#### launching R calculation ########"
Rscript /home/lescailab/CODE/hgtseq_tools/calculate_heaptmap.R $datafile $parallel $SLURM_CPUS_PER_TASK