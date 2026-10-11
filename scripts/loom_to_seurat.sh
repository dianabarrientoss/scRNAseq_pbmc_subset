#!/bin/bash
#SBATCH --job-name=loom_to_seurat
#SBATCH --partition=defq
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/out_logs/loom_to_seurat_%j.log

module load apptainer/1.5.3

echo "Nodo: $(hostname)"
echo "Inicio: $(date)"

USER_DATA_DIR="/mnt/data/bioinfo3/$USER"
TMPDIR=$(mktemp -d)
mkdir -p $TMPDIR/tmp

cd $USER_DATA_DIR/proyecto_integrador/scripts

apptainer exec \
  --home "$USER_DATA_DIR:/home/$USER" \
  -B $USER_DATA_DIR \
  -B $TMPDIR/tmp:/tmp \
  -B /cm/shared/apps/r/4.6.0-studio/lib64/R/library:/opt/curso_lib \
  -B /cm/local/apps/gcc/14.2.0/lib64:/opt/gcc14/lib64 \
  -B /cm/shared/apps/r/4.6.0-studio/lib64/R:/opt/R_base \
  -B /cm/shared/apps/hdf5/1.14.0/lib:/opt/hdf5/lib \
  --env LD_LIBRARY_PATH="/opt/hdf5/lib:/opt/R_base/lib:/opt/gcc14/lib64:$LD_LIBRARY_PATH" \
  --env R_LIBS_SITE="/opt/curso_lib" \
  --env LANG=en_US.UTF-8 \
  --env LC_ALL=en_US.UTF-8 \
  /cm/shared/apps/rstudio/rstudio_sandbox \
  Rscript loom_to_seurat_V3.r

echo "Fin: $(date)"
