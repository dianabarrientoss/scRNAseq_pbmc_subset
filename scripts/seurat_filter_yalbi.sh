#!/bin/bash
#SBATCH --job-name=seurat_filter
#SBATCH --partition=defq
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/out_logs/seurat_filter_%j.log
#SBATCH --error=/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/out_logs/seurat_filter_%j.err

module load apptainer/1.5.3
module load r/4.6.0

Rscript /mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_filter_yalbi.R
