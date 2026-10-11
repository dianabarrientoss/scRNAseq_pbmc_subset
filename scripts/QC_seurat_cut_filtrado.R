
# DE FORMATO LOOM A FORMATO SEURAT
# AUTORES:  "Barrientos-González Diana, Cervantes-Levario Román & Contreras-Rivera Jair"
# FECHA: "23-09-2026"

setwd("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts")

seurat_cut_filtrado <- readRDS("seurat_cut_filtrado.rds")

dim(seurat_cut_filtrado)                     # genes x células
Assays(seurat_cut_filtrado)                  # RNA, SCT, etc.
Layers(seurat_cut_filtrado[["RNA"]])         # counts, data, scale.data
Reductions(seurat_cut_filtrado)              # pca, umap, tsne si existen
colnames(seurat_cut_filtrado@meta.data)      # todas las columnas de metadatos
str(seurat_cut_filtrado, max.level = 2)      # estructura completa
summary(seurat_cut_filtrado$nCount_RNA)      # distribución de una métrica
table(seurat_cut_filtrado$orig.ident)        # células por muestra
head(rownames(seurat_cut_filtrado))          # nombres de genes




seurat_cut_filtrado[["percent.mt"]] <- PercentageFeatureSet(seurat_cut_filtrado, pattern = "^MT-")

VlnPlot(seurat_cut_filtrado, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, col = "pink")


