library(Seurat)
library(DropletUtils)

ruta_data_yalbi <- file.path("/mnt/data/bioinfo3", "dbarrientos", "proyecto_integrador", "scripts", "proyecto-equipo-B.rds")
seurat_yalbi <- readRDS(ruta_data_yalbi)
seurat_yalbi <- as.Seurat(seurat_yalbi, counts = "counts", data = NULL)
saveRDS(seurat_yalbi, file = "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_yalbi.rds")

