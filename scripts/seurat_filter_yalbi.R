library(Seurat)
library(Matrix)

seurat_obj <- readRDS("/home/jcontreras/proyecto/scripts/seurat_yalbi.rds")
counts_all <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

tot <- Matrix::colSums(counts_all)

cat("Diagnostico\n")
cat("Barcodes totales:", length(tot), "\n")
cat("Barcodes con 0 UMI:", sum(tot == 0), "\n")
cat("Barcodes con >0 UMI:", sum(tot > 0), "\n")
cat("UMI maximo:", max(tot), "\n")
cat("UMI minimo no-cero:", min(tot[tot > 0]), "\n")
print(quantile(tot[tot > 0], c(0, .5, .9, .95, .99, 1)))

cutoff <- 0.10 * max(tot)
cat("\nCutoff 10% maximo =", round(cutoff, 1), "UMI\n")
cat("Barcodes >= cutoff:", sum(tot >= cutoff),
    sprintf("(%.2f%%)\n", 100 * sum(tot >= cutoff) / length(tot)))
cat("Barcodes < cutoff:", sum(tot < cutoff),
    sprintf("(%.2f%%)\n", 100 * sum(tot < cutoff) / length(tot)))

out_dir <- "/home/jcontreras/proyecto/results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf_cutoff <- file.path(out_dir, "barcode_rank_yalbi_max10.pdf")

pdf(pdf_cutoff, width = 9, height = 6)
tot_sorted <- sort(tot, decreasing = TRUE)
rank <- seq_along(tot_sorted)
nz <- tot_sorted > 0

plot(rank[nz], tot_sorted[nz], log = "xy", type = "l", lwd = 2,
     xlab = "Rank de barcode", ylab = "UMI totales (log)",
     main = "Barcode rank yalbi - cutoff 10% maximo")
abline(h = cutoff, col = "red", lty = 2, lwd = 2.5)
abline(h = c(100, 500, 1000), col = "gray60", lty = 4, lwd = 1)
legend("bottomleft",
       legend = c(sprintf("Cutoff 10%% max = %.0f UMI", cutoff),
                  "100 UMI", "500 UMI", "1000 UMI"),
       col = c("red", "gray60", "gray60", "gray60"),
       lty = c(2, 4, 4, 4), lwd = c(2.5, 1, 1, 1),
       bg = "white", cex = 0.9)
dev.off()

cat("PDF guardado:", pdf_cutoff,
    "(", round(file.info(pdf_cutoff)$size / 1024, 1), "KB )\n")

keep <- tot >= cutoff
counts_filtrado <- counts_all[, keep]

cat("\nFiltrado\n")
cat("Celulas antes:", ncol(counts_all), "\n")
cat("Celulas despues:", ncol(counts_filtrado), "\n")
cat("Retenidas:", sprintf("%.2f%%\n",
    100 * ncol(counts_filtrado) / ncol(counts_all)))

pdf_ba <- file.path(out_dir, "umi_before_after_yalbi.pdf")
pdf(pdf_ba, width = 10, height = 5)
par(mfrow = c(1, 2))

tot_pos <- tot[tot > 0]
hist(log10(tot_pos + 1), breaks = 60, col = "gray70", border = "white",
     main = sprintf("Antes (n=%d con UMI)", length(tot_pos)),
     xlab = "log10(UMI + 1)", ylab = "Frecuencia")

tot_after <- Matrix::colSums(counts_filtrado)
if (length(tot_after) > 0) {
    hist(log10(tot_after + 1), breaks = 60, col = "forestgreen",
         border = "white",
         main = sprintf("Despues (n=%d)", length(tot_after)),
         xlab = "log10(UMI + 1)", ylab = "Frecuencia")
} else {
    plot.new()
    text(0.5, 0.5, "Sin celulas tras filtrado", cex = 1.3, col = "red")
}
par(mfrow = c(1, 1))
dev.off()

cat("PDF guardado:", pdf_ba,
    "(", round(file.info(pdf_ba)$size / 1024, 1), "KB )\n")

out_data <- "/home/jcontreras/proyecto/scripts"
dir.create(out_data, showWarnings = FALSE, recursive = TRUE)

saveRDS(counts_filtrado, file.path(out_data, "counts_filtrado_yalbi.rds"))

seurat_filtrado <- CreateSeuratObject(counts = counts_filtrado,
                                      project = "yalbi_filtrado")
saveRDS(seurat_filtrado,
        file.path(out_data, "seurat_yalbi_filtrado.rds"))

