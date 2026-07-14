#!/bin/sh
#SBATCH --partition=HMEM
#SBATCH --nodes=1
#SBATCH --mail-user=aabdala@nioz.nl
#SBATCH --mail-type=FAIL,BEGIN,END,TIME_LIMIT_80
#SBATCH --time=5-10:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=528G    #Default Mem

module load anaconda/2024.02
conda activate /export/lv10/projects/projects_WR/envs/Meta-cascabel-test

# conda activate /export/lv1/user/aabdala/.conda/envs/metac4/

snakemake --configfile config.yaml  -j3 -c100 --keep-going 
# --use-conda --conda-frontend conda
# --rerun-triggers mtime
snakemake --configfile config.yaml --report report.hpc.zip

