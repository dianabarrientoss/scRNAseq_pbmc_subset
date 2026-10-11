# DE FORMATO LOOM A FORMATO SEURAT
# AUTORES:  "Barrientos-González Diana, Cervantes-Levario Román & Contreras-Rivera Jair"
# FECHA: "23-09-2026"

# Objetivo: crear un loom pequeño y estándar para poder probarlo


# Cores disponibles
n_cores <- parallel::detectCores()

# Dejar 2 cores libres para el sistema
options(mc.cores = max(1, n_cores - 2))

cat("Cores detectados:", n_cores, "\n")
cat("Cores configurados para R:", getOption("mc.cores"), "\n")
setwd("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts")

file <- "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/proyecto-equipo-B.loom"

library(loomR)

# Abrir archivo original:
lfile <- loomR::connect(
  filename = file,
  mode = "r",
  skip.validate = TRUE
)

# Número de células para prueba:
n_use <- 1000
cells_use <- 1:n_use

# Extraer solo esas células
# loomR presenta matrix como células x genes
mat_sub <- lfile[["matrix"]][cells_use, ]

# Metadatos de genes (armar una lista con ID, Symbol, Type y rownames de los genes):
gene_attrs <- list(
  ID       = as.character(lfile$row.attrs$ID[]),
  Symbol   = as.character(lfile$row.attrs$Symbol[]),
  Type     = as.character(lfile$row.attrs$Type[]),
  rownames = as.character(lfile$row.attrs$rownames[])
)

#Metadatos de células
#lista con los barcodes de las 1,000 células elegidas.
cell_attrs <- list(
  colnames = as.character(
    lfile$col.attrs$colnames[cells_use]
  )
)

# Crear un loom nuevo y estándar
small_loom <- loomR::create(
  filename = "proyecto-equipo-B_1k.loom",
  data = mat_sub,
  gene.attrs = gene_attrs,
  cell.attrs = cell_attrs,
  do.transpose = FALSE,
  overwrite = TRUE
)

small_loom$close_all()
lfile$close_all()
