library(BiocManager)

.libPaths(c(
    "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    .libPaths()
))

BiocManager::install(
    c(
        "SingleCellExperiment",
        "Matrix",
        "DropletUtils",
        "SoupX"
    ),
    lib = "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    ask = FALSE,
    update = FALSE
)
