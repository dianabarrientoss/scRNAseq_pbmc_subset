#!/bin/bash

#SBATCH --job-name=seurat_verificacion
#SBATCH --partition=defq
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=../out_logs/loom_%j.out
#SBATCH --error=../out_logs/loom_%j.err

module load r/4.4.1-ybalderas

Rscript --vanilla - <<'EOF'
.libPaths(c(
    "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    .libPaths()
))

library(Seurat)
packageVersion("Seurat")
