library(BiocManager)

.libPaths(c(
    "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    .libPaths()
))

BiocManager::install(
    "rhdf5",
    lib = "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    ask = FALSE,
    update = FALSE
)

BiocManager::install(
    "LoomExperiment",
    lib = "/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/libraries/Rlibs",
    ask = FALSE,
    update = FALSE
)
