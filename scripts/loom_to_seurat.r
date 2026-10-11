# DE FORMATO LOOM A FORMATO SEURAT
# AUTORES:  "Barrientos-González Diana, Cervantes-Levario Román & Contreras-Rivera Jair"
# FECHA: "23-09-2026"

setwd("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts")

library(rhdf5)
library(SingleCellExperiment)
library(Seurat)
library(SeuratDisk)

#Versiones con las que arranca el job
print(sessionInfo())

#Conectar al loom
ruta_loom <- "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/data/proyecto-equipo-B.loom"
mi_loom <- Connect(filename = ruta_loom, mode = 'r')

#Nombres de genes y células
sym     <- mi_loom[["row_attrs/Symbol"]][]  #spimbolos de genes (sus nombres)
rn      <- mi_loom[["row_attrs/rownames"]][]  #IDs para ensembl
celulas <- mi_loom[["col_attrs/colnames"]][]  #barcodes de células

#Symbol tiene duplicados que seurat no acepta, pero los de ensembl no tienen duplicados. Por eso ensembl es el identificador de filas.
cat("Símbolos duplicados:", sum(duplicated(sym)), "\n")
cat("Ensembl duplicados:", sum(duplicated(rn)), "\n")

n_gen <- length(rn)
n_cel <- length(celulas)
cat("Genes:", n_gen, " Células:", n_cel, "\n")

#Leer la matriz por bloques y convertir a dispersa

mat  <- mi_loom[["matrix"]]
dims <- mat$dims
cat("Dimensiones de /matrix:", dims, "\n")

#esto permite comprobar que coinciden con el número de genes y células
stopifnot(setequal(dims, c(n_cel, n_gen))) 
celulas_en_filas <- (dims[1] == n_cel)

#leer por bloques evita cargar todos los datos de una
tam_bloque <- 1000
bloques <- split(seq_len(n_cel), ceiling(seq_len(n_cel) / tam_bloque))
lista <- vector("list", length(bloques))
t0 <- Sys.time()
#Leer cada bloque, transponer si hace falta y convertirlo a matriz dispersa
for (i in seq_along(bloques)) {
  idx <- bloques[[i]]
  m <- if (celulas_en_filas) t(mat[idx, ]) else mat[, idx]
  lista[[i]] <- Matrix::Matrix(m, sparse = TRUE) #es la de genes por células
  rm(m)
  
  if (i == 1) {
    nnz_est <- length(lista[[1]]@x) * length(bloques)
    cat("No ceros en el bloque 1:", length(lista[[1]]@x),
        "| estimado total:", format(nnz_est, big.mark = ","), "\n")
    # Una matriz dispersa de R admite alrededor de 2.147 millones de no ceros
    if (nnz_est > 2.1e9) {
      stop("El estimado de no ceros excede el límite de la matriz dispersa. ",
           "Hay que hacer subconjunto de células.")
    }
  }
  if (i %% 10 == 0 || i == 1) {
    cat("Bloque", i, "de", length(bloques), "| minutos:",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
  }
}
#unir los bloques por columnas (células). 
#después: se eliminan objetos intermedios, liberar memoria
matriz_cuentas <- do.call(cbind, lista)
rm(lista); gc()

cat("Clase:", class(matriz_cuentas), " Dim:", dim(matriz_cuentas), "\n")

#asignar nombre de columnas y filas
rownames(matriz_cuentas) <- rn
colnames(matriz_cuentas) <- celulas

#crear tabla de equivalencias ensembl: símbolo original y símbolo único
tabla_genes <- data.frame(
  ensembl      = rn,
  symbol       = sym,
  symbol_unico = make.unique(sym),
  row.names    = rn
)

#Cerrar el loom cuando ya se leyó todo lo que necesitamos:
mi_loom$close_all()

#Guardar la tabla de genes por si algo falla después
dir.create("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/output",
           showWarnings = FALSE, recursive = TRUE)

saveRDS(
  tabla_genes,
  file = "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/output/tabla_genes.rds"
)

#crear objeto Seurat a partir de la matriz de conteos sin normalizar

mi_seurat <- CreateSeuratObject(counts = matriz_cuentas, project = "proyecto-equipo-B")
rm(matriz_cuentas); gc()

#intentar añadir los símbolos de genes como anotación
tryCatch(
  mi_seurat[["RNA"]][["symbol"]] <- tabla_genes[rownames(mi_seurat), "symbol_unico"],
  error = function(e) message("No se pudo agregar 'symbol' al assay: ", conditionMessage(e))
)

#Revisión de las dimensiones, las anotaciones
print(mi_seurat)
print(dim(mi_seurat))
print(head(tabla_genes))
print(GetAssayData(mi_seurat, layer = "counts")[1:5, 1:5])

#Guardar el objeto
saveRDS(
  mi_seurat,
  file = "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/output/proyecto-equipo-B_seurat.rds"
)

cat("Listo:", format(Sys.time()), "\n")

#Información de la sesión
print(sessionInfo())

writeLines(capture.output(sessionInfo()),
           "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/output/sessionInfo_loom_to_seurat.txt")

#se utilizó como referencia:
# https://github.com/satijalab/seurat/blob/master/vignettes/conversion_vignette.Rmd
# lineas 104-110 de código pero no todas se usan
