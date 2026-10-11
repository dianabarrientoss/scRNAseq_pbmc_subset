############ Configuraciones iniciales ############################

# Semilla
set.seed(20260827)

# Libreria
library(Seurat)
library(DropletUtils)

# Usuario
#usuario <- Sys.getenv("USER")
usuario <- "dbarrientos"

###################################################################


############### Funcion para contestar las preguntas ##############

adress_answers <- function(seurat, file_type){
	cat("\n El documento es el:", file_type, "\n")
	
	cat("\n El tipo de objeto para", file_type, "es", class(seurat)[1], "\n")
	
	cat("\n El número de genes x células (primero genes, segundo células): ", nrow(seurat), "y", ncol(seurat), "\n")
	
	cat("\n Las capas del assay de RNA son:", paste(Layers(seurat[["RNA"]]), collapse = ", "), "\n")
	
	cat("\n Los genes tienen los atributos:", paste(colnames(seurat[["RNA"]][[]]), collapse = ", "), "\n")
	
	# Vamos a tomar una submuestra para evitar calcular todas las estadísticas
	
	idx <- sample(ncol(seurat), min(5000, ncol(seurat)))
	
	sub <- as.matrix(GetAssayData(seurat, assay = "RNA", layer = "counts")[, idx])
	
	if (all(sub >= 0)) {cat("Todos los conteos son no negativos\n")} 
	else {cat("Existen conteos negativos\n")}
	
	if (all(sub == round(sub))) {cat("Todos los conteos son enteros\n")} 
	else {cat("Existen conteos decimales\n")}
	
	cat("Summary del tamaño por célula:\n")
	
	print(summary(Matrix::colSums(sub)))

	assay_use <- DefaultAssay(seurat)	

	counts_all <- GetAssayData(seurat, assay = assay_use, layer = "counts")
	
	tot  <- Matrix::colSums(counts_all)

	tot  <- sort(tot, decreasing = TRUE)
	
	rank <- seq_along(tot)

	cat("N barcodes:", length(tot),"| zeros:", sum(tot == 0),"| min:", min(tot), "| max:", max(tot), "\n")

	keep <- tot > 0 & !is.na(tot)

	out_pdf <- file.path("/mnt/data/bioinfo3", usuario, "proyecto_integrador", "results", "figures", paste0("barcode_ranks_", file_type, ".pdf"))

	pdf(out_pdf, width = 7, height = 6)
	
	plot(rank[keep], tot[keep], log = "xy", type = "l", lwd = 2, xlab = "Rank de barcode", ylab = "UMI totales", main = paste("Barcode rank plot —", file_type))
	
	abline(h = c(100, 500, 1000), col = c("orange","blue","red"), lty = 2)

	legend("bottomleft", legend = c("100","500","1000"), col = c("orange","blue","red"), lty = 2, title = "UMI umbral")
	
	dev.off()
	
	# Metadatos 
	cd <- seurat@meta.data

	cand <- grep("batch|lote|donor|donante|patient|subject|sample|muestra|condition|cond|treat|group", names(cd), ignore.case = TRUE, value = TRUE)

	cat("\nColumnas de metadatos candidatas:\n")
	print(cand)

	cat("\nNúmero total de células:", ncol(seurat), "\n")
	
	cat("Número total de genes :", nrow(seurat), "\n")
	
}


ruta_data_cut <- file.path("/mnt/data/bioinfo3", usuario, "proyecto_integrador", "scripts", "proyecto-equipo-B_10000cells_seurat.rds")
seurat_cut <- readRDS(ruta_data_cut)
adress_answers(seurat_cut, "File cut")
