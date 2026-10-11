# DE FORMATO LOOM A FORMATO SEURAT
# AUTORES:  "Barrientos-González Diana, Cervantes-Levario Román & Contreras-Rivera Jair"
# FECHA: "23-09-2026"


# instalar scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# instalar loomR de GitHub usando: remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')library(loomR)
library(Seurat)
library(SeuratDisk)
library(patchwork)

# Usar devtools para instalar hdf5r y loomR desde Github


# Instalar Devtools desde CRAN
install.packages("devtools")
# Usar devtools para instalar hdf5r y loomR desde github:
devtools::install_github(repo = "hhoeflin/hdf5r")

devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

library(loomR)


setwd("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts")

file <- "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/proyecto-equipo-B.loom"

pbmc_loom <- read(file)

LoadLoom(
  file = file,
  
)

library(loomR)

#abrir loom saltando la validación estricta
lfile <- loomR::connect(
  filename = "proyecto-equipo-B.loom",
  mode = "r",
  skip.validate = TRUE
)

lfile


names(lfile$row.attrs)
names(lfile$col.attrs)

# Extraer nombres de features y células
features <- lfile$row.attrs$rownames[]
cells <- lfile$col.attrs$colnames[]

# Revisar tipos de features antes de crear Seurat
print(table(lfile$row.attrs$Type[]))

# Extraer matriz del loom
count_matrix <- lfile$matrix[, ]

# loomR la presenta como células x features,
# Seurat necesita features x células así que se transpone la matriz t()
count_matrix <- t(count_matrix)

# Convertir a matriz sparse
#esto permite eliminar los datos que son cero y que sea más manejable la matriz
counts <- Matrix(counts, sparse = TRUE)

# Asignar nombres
# make.unique() evita nombres duplicados heciendo identificadores únicos
rownames(counts) <- make.unique(as.character(features))
#asignar estos nuevos nombres a las columnas
colnames(counts) <- make.unique(as.character(cells))

# Crear objeto Seurat
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "proyecto-equipo-B"
)

# Mostrar objeto para comprobar que se creó correctamente
print(seurat_obj)


#checar columnas y nombres. Esto sí funcionó para checar qué contiene y cómo sse llama todo
library(hdf5r)

f <- H5File$new(file, mode = "r")

# Ver qué hay en row_attrs (atributos de genes)
print(f[["row_attrs"]]$names)

# Ver qué hay en col_attrs (atributos de células)
print(f[["col_attrs"]]$names)

f$close_all()


library(Seurat)
library(loomR)

# Abrir el loom, 
lfile <- loomR::connect(
  filename = file,
  mode = "r",
  skip.validate = TRUE
)

# Convertir directamente a Seurat (no funcionó)
seurat_obj <- as.Seurat(
  lfile, #llamo a LOOM de arribita
  assay = NULL,
  cells = "col_names",
  features = "ID",
)

# Revisar
seurat_obj

#Ahora crear subset para que no sea el atlas de células:

library(loomR)

loom <- connect(file, mode = "r", skip.validate = TRUE)

# Leer solo las primeras 1000 células
mat_sub <- loom[["matrix"]][, 1:1000]

# Genes
genes <- loom$row.attrs$col_names[]


loom$close_all()


#Checar la memoira:

# Ver RAM total del sistema
system("sysctl -n hw.memsize") |> as.numeric() / 1024^3

# Aumentar el tamaño máximo de objetos globales
options(future.globals.maxSize = 16 * 1024^3)  # 16 GB, ajustalo según tu RAM

# Para Seurat específicamente
options(Seurat.object.assay.version = "v5")


# Poner esto al inicio de tu script, antes de cargar el archivo
#Esto ayuda a evitar errores cuando Seurat necesita trabajar con objetos grandes
options(future.globals.maxSize = 16 * 1024^3)



#SI SOLUCIONÓ PROBLEM  GB DE MATRIX
seurat_obj <- SeuratDisk::LoadLoom(
  file = file,
  assay = "RNA",
  cells = "colnames",
  features = "Symbol",
  filter = "none",
  verbose = TRUE
)

seurat_obj

#guardar objeto seurat como archivo loomn:
pbmc.loom <- as.loom(pbmc, filename = ruta_loom, verbose = FALSE)
pbmc.loom


#Convertir nuevamente el objeto loom a seurat:
pbmc.seurat <- as.Seurat(pbmc.loom)


#visualizar con violin plot la distribución de expresión de cuatro genes
#Sparc, Ftl1, Junb y Ccl4 se utilizaron como ejemplo  para inspeccionar sus patrones de expresión en las células de este dataset

VlnPlot(pbmc.seurat, features = c("Sparc", "Ftl1", "Junb", "Ccl4"), ncol = 2, pt.size = 0.1)
