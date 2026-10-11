# DE FORMATO LOOM A FORMATO SEURAT
# AUTORES:  "Barrientos-González Diana, Cervantes-Levario Román & Contreras-Rivera Jair"
# FECHA: "23-09-2026"


#Configuración:

# Cores visibles para R
n_cores <- parallel::detectCores()

# Dejar 2 libres
options(
  mc.cores = max(1, n_cores - 2)
)

cat("Cores detectados:", n_cores, "\n")
cat("Cores configurados:", getOption("mc.cores"), "\n")

setwd("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts")

file <- "proyecto-equipo-B.loom"

library(loomR)


#Abrir loom originaL:

lfile <- loomR::connect(
  filename = file,
  mode = "r",
  skip.validate = TRUE
)


#Se crea subset de las primeras mil células:

n_use <- 1000L

# Seguridad: es importante no pedir más células de las que existen:
# Determinar el número total de células disponibles en el archivo Loom  y seleccionar cuántas se utilizarán para el análisis.
n_total <- length(lfile$col.attrs$colnames[])

n_use <- min(n_use, n_total)

#Generar los índices de las primeras n_use células que serán extraídas
cells_use <- seq_len(n_use)

#número total de células y número seleccionado
cat("Células totales:", n_total, "\n")
cat("Células a extraer:", n_use, "\n")


#Leer la matriz, imprimir dim() y su tamaño en RAM:


cat("Leyendo matriz...\n")

mat_sub <- lfile[["matrix"]][cells_use, ]

cat("Dimensión extraída:\n")
print(dim(mat_sub))

cat(
  "Tamaño en RAM:",
  format(object.size(mat_sub), units = "auto"),
  "\n"
)


#Metadatos mínimos necesarios:
#Se conserva solo symbol para genes y colnames para células
# se quitaron ID, type y rownames para descartar que el problema estuviera en los metadatos
gene_attrs <- list(
  Symbol = as.character(
    lfile$row.attrs$Symbol[]
  )
)

cell_attrs <- list(
  colnames = as.character(
    lfile$col.attrs$colnames[cells_use]
  )
)


#Creear loom reducido (El nombre sale de n_use)

out_file <- paste0(
  "proyecto-equipo-B_",
  n_use,
  "cells.loom"
)

cat("Creando:", out_file, "\n")

small_loom <- loomR::create(
  filename = out_file,
  data = mat_sub,
  gene.attrs = gene_attrs,
  cell.attrs = cell_attrs,
  do.transpose = FALSE,
  overwrite = TRUE,
  display.progress = TRUE
)


#Cerrar archivos y liberar memoria:

small_loom$close_all()
lfile$close_all()

rm(
  small_loom,
  lfile,
  mat_sub,
  gene_attrs,
  cell_attrs
)

gc()

cat("\nListo:", out_file, "\n")
