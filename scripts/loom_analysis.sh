#!/bin/bash

#SBATCH --job-name=loom_install
#SBATCH --partition=defq
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=../out_logs/loom_analysis_%j.out
#SBATCH --error=../out_logs/loom_analysis_%j.err

module load r/4.4.1-ybalderas
cd /mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/

Rscript loom_analysis.R

