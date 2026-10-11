#!/bin/bash

#SBATCH --job-name=seurat-disk_install
#SBATCH --partition=defq
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=../../out_logs/loom_%j.out
#SBATCH --error=../../out_logs/loom_%j.err

module load r/4.4.1-ybalderas
moduel load hdf5/1.14.0

Rscript --vanilla - <<'EOF'
.libPaths(c(
    "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    .libPaths()
))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("mojaveazure/seurat-disk")
EOF

#mojaveazure.github.io/seurat-disk en README


