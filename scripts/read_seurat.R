
#Este archivo permite leer el .rds que preparó la profesora y obtener un objeto Seurat


#Leer el seurat:
data <- readRDS("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/proyecto-equipo-B.rds")

#cargar la librería necesaria:
library(Seurat)

#convierte el objeto de entrada en seurat. toma la matriz de conteos crudos.
data_seurat <- as.Seurat(
  data,
  counts = "counts",
  data = NULL      #esto indica que non hay matriz normalizada
)

#imprimir el resumen:
data_seurat
