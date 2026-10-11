#Fijamos semilla
set.seed(20260827)


## Carga de librerás para los archivos
## Definimos la ruta a donde descargamos las libreriás
.libPaths(c(
    "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    .libPaths()
))


#Y las cargamos
library(LoomExperiment)
library(rhdf5)
library(DropletUtils)


#Cargamos la ruta de nuestro archivo loom
ruta_data <- "/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-B.loom"


#Loo leemos y exploramos los artibutos
loom <- LoomExperiment::import(ruta_data)
loom


#Para ver el tipo de nuestro archivo, si es seurat o single cell experimento ocupamos leer los atributos con la librería rhdf5
global_attrs <- h5readAttributes(ruta_data, name = "/")
cat("\nEl tipo de experimento usado es:", global_attrs$`LoomExperiment-class`, "\n")


#Para ver el numero de genes x células
cat("\nEl número de genes x células (El primer valor corresponde al número de genes y el segundo al de células): ", dim(loom), "\n")


#Para explorar el nombre de la matriz de conteos:
cat("\n La matriz de conteos se llama:", assayNames(loom), "\n")


#Para ver los atributos de los genes
cat("\n Los genes tienen los atributos:", names(rowData(loom)),  "\n")


#Para revisar si las cuentas vienen crudas o no, dado que es costoso computacionalmente tomamos una submuestra de 5,000
idx <- sample(ncol(loom), min(5000, ncol(loom)))
sub <- as.matrix(assay(loom, "counts")[, idx])

if (all(sub >= 0)) cat("Todos los conteos son no negativos\n") else cat("Existen conteos negativos\n")

if (all(sub == round(sub))) cat("Todos los conteos son enteros\n") else cat("Existen conteos decimales\n")

cat("Resumen del tamaño de librería por célula:\n")
print(summary(colSums(sub)))


# Para revisar si las cuentas están filtradas o no intentarmeos filtrar con dropletutilities, y verificaremos la curva
br <- barcodeRanks(assay(loom, "counts"))
knee <- metadata(br)$knee
inflection <- metadata(br)$inflection

pdf("/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/results/figures/barcode_ranks.pdf", width = 7, height = 6)

plot(
    br$rank,
    br$total,
    log = "xy",
    pch = 16,
    cex = 0.3,
    xlab = "Rank del barcode",
    ylab = "UMIs totales",
    main = "Barcode rank plot"
)

o <- order(br$rank)

lines(
    br$rank[o],
    br$fitted[o]
)

abline(h = knee, lty = 2)
abline(h = inflection, lty = 2)

legend(
    "bottomleft",
    legend = c("Knee", "Inflection"),
    lty = 2,
    bty = "n"
)

dev.off()


# Para revisar si hay informacion importante de lote, donante y condicion
cd <- as.data.frame(colData(loom))

cand <- grep(
    "batch|lote|donor|donante|patient|subject|sample|muestra|condition|cond|treat|group",
    names(cd),
    ignore.case = TRUE,
    value = TRUE
)

print(cand)
print(nrow(loom))
