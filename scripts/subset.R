library(hdf5r)
library(Matrix)
library(Seurat)

#definir el archivo y cuántas células se quieren
file <- "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/proyecto-equipo-B.loom"

n_use <- 10000L

#abrir el loom 
h5 <- hdf5r::H5File$new(
  file,
  mode = "r"
)

# Dataset de expresión
mat_h5 <- h5[["matrix"]]

# importante:
# En hdf5r usar $dims, NO dim()
#se guarda en "d" las dimensiones del dataset matrix
d <- mat_h5$dims

cat(
  "Dimensiones de /matrix:",
  paste(d, collapse = " x "),
  "\n"
)

#leer los símbolos de genes y los barcodes (identificador de células) para contar cuántos hay:

genes <- as.character(
  h5[["row_attrs"]][["Symbol"]][]
)

cells <- as.character(
  h5[["col_attrs"]][["colnames"]][]
)

cat("Features:", length(genes), "\n")
cat("Células:", length(cells), "\n")

#seleccionar células:

#para evitar pedir más células de las que existen, toma el valor de 10000L que fijamos al inicio
n_use <- min(
  n_use,
  length(cells)
)

#tomar las primeras en el orden del archivo:
cells_use <- seq_len(n_use)


cat("Células a cargar:", n_use, "\n")

#Leer solo el subset:
#compara "d" con el número de genes y células. Si la matriz es genes × células, toma las columnas. Si es células x genes, lo transpone --> t(counts): fue nuestro caso
#la última parte ocurre cuando no coincide con ninguno de los dos y se marca como error

if (
  d[1] == length(genes) &&
  d[2] == length(cells)
) {
  
  cat("Orientación: features x células\n")
  
  counts <- mat_h5[
    ,
    cells_use
  ]
  
} else if (
  d[1] == length(cells) &&
  d[2] == length(genes)
) {
  
  cat("Orientación: células x features\n")
  
  counts <- mat_h5[
    cells_use,
  ]
  
  counts <- t(counts)
  
} else {
  
  stop(
    paste0(
      "Dimensiones inesperadas: ",
      paste(d, collapse = " x "),
      "\nFeatures: ", length(genes),
      "\nCélulas: ", length(cells)
    )
  )
}


#cerrar HDF5:

h5$close_all()

#Convertir a sparse:
#esto guardará solo los conteos que sean distintos de cero para hacer más manejable la matriz
cat("Convirtiendo a sparse...\n")

counts <- Matrix::Matrix(
  counts,
  sparse = TRUE
)

#Asignar nombres únicos:
#los sínbolos repetidos de genes en las filas se pondrían como: gen, gen.1, gen.2, etc.
#es solo una etiqueta técnica
rownames(counts) <- make.unique(genes)

colnames(counts) <- make.unique(
  cells[cells_use]
)

#imprimir los 38606 features × 10000 células y el tamaño en RAM:
cat(
  "Matriz final:",
  nrow(counts),
  "features x",
  ncol(counts),
  "células\n"
)

cat(
  "Tamaño en RAM:",
  format(object.size(counts), units = "auto"),
  "\n"
)

# Crear seurat:
#se crea con "proyecto-equipo-B" como proyecto, mínimo de células y features se deja como cero.
#Se crean automáticamente los metadatos nCount_RNA y nFeature_RNA

seurat_obj <- Seurat::CreateSeuratObject(
  counts = counts,
  project = "proyecto-equipo-B",
  min.cells = 0,
  min.features = 0
)

print(seurat_obj)

#Guardar el .rds:
#guarda el archivo en en el directorio de trabajo de la sesión

output_file <- paste0(
  "proyecto-equipo-B_",
  n_use,
  "cells_seurat.rds"
)

saveRDS(
  seurat_obj,
  output_file
)

cat("\nGuardado en:", output_file, "\n")


# Ver QC metrics en VlnPlot
#Esto permite ver la distribución
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, col ="pink")

