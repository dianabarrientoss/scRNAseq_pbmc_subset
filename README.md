# Proyecto final — equipo B
> Plantilla del módulo de análisis de células individuales, LCG 2027-1.
> Equipo: Diana Barrientos González, Román Cervantes Levario, Jair Emiliano Contreras Rivera
## Qué es este repositorio
El análisis del conjunto de datos que se le asignó a tu equipo, de principio a
fin, de forma que **cualquier persona pueda reproducirlo desde cero**. 
Por su tamaño, no se pueden cargar los datos completos (927,681 células × 38,606 números). Por eso hay varios scripts: unos intentan la conversión completa, otros trabajan con un subconjunto de células, y uno lee el objeto ya recortado
## Integrantes
| Nombre | Qué hizo |
|---|---|
|Diana Barrientos| Conversion de archivo loom a seurat, Genercación de subset de datos, Quality control posteior al filtrado, Correccion de issues presentados, PCA, anotacion y UMAP, recreacion de resultados |
|Roman Cervantes | Readme, Redaccion de preguntas, Estructura del proyecto de github, Interpretacion biologica de datos, Reporte e instrucciones de como correr los datos, Interpretacion de figuras |
|Jair Rivera |  Instalación de librerias, Exploracion y analisis de las preguntas iniciales, filtrado de datos, correccion de errores relacionados al filtrado, reporte de recreacion de resultados |

## Contexto biologico


Los datos corresponden a PBMCs aisladas y criopreservadas por AllCells. Las células son primarias y contienen aproximadamente 1 pg de ARN por célula. Para la preparación de las bibliotecas se utilizó el protocolo Chromium Single Cell 3′ Reagent Kits v3.1 de 10x Genomics. Se prepararon ocho bibliotecas en un mismo chip, del cual se analizó el canal 5. Se esperaban 5,000 células, pero se detectaron 4,053. La secuenciación se realizó en Illumina NovaSeq, con aproximadamente 52,546 lecturas por célula. Las lecturas incluyeron los códigos de barras celulares y los UMI en Read 1, mientras que Read 2 contenía las secuencias de los transcritos. Para el procesamiento se utilizó el parámetro `--expect-cells=5000`. Los datos cuentan con licencia CC BY 4.0.

## Cómo reproducir este análisis

### Paso 0. Abrir RStudio Server en el clúster (desde el navegador)
Todo el código se corre dentro del clúster, usando un RStudio Server que se lanza con un script y se abre en el navegador de tu laptop. Necesitas dos terminales abiertas al mismo tiempo.
> IMPORTANTE: copia el script en tu directorio. Cada persona debe copiarlo a su propio directorio y correrlo desde su copia.
> El script original está en la carpeta de Diana:
> `/mnt/data/bioinfo3/dbarrientos/lanzar_rstudiov3.sh’

```bash
# 1) Copiar el script a tu directorio (reemplazando “tu_usuario”)
cp /mnt/data/bioinfo3/dbarrientos/lanzar_rstudiov3.sh /mnt/data/bioinfo3/tu_usuario/
```

Nota: aunque dice que hay que leer los datos sin copiarlos, en nuestro caso sí fue necesario copiar los datos.
#### Terminal 1 — lanzar el servidor (dentro del clúster)

Conéctate al clúster, crea una sesión interactiva, abre `bash` y ejecuta tu copia del script:

```bash
bash
bash /mnt/data/bioinfo3/TU_USUARIO/lanzar_rstudiov3.sh
```

La terminal mostrará algo como esto (los valores cambian según tu usuario y el nodo asignado):

```text
¡RStudio Server corriendo en el nodo: node01!
Directorio de trabajo: /mnt/data/bioinfo3/TU_USUARIO
Puerto local del clúster: 8794
----------------------------------------------------
```
También, es normal que aparezcan mensajes como:

TTY detected. Printing informational message about logging configuration. Logging configuration loaded from '/etc/rstudio/logging.conf'. Logging to '/home/dbarrientos/.local/share/rstudio/log/rserver.log'.
TTY detected. Printing informational message about logging configuration. Logging configuration loaded from '/etc/rstudio/logging.conf'. Logging to '/home/dbarrientos/.local/share/rstudio/log/rsession-dbarrientos.log'.

> Es importante que no corran nada en esta terminal ni en la segunda (a continuación), debes dejarlo como está y abrir una tercera ventana para trabajar.

#### Terminal 2 — (en tu laptop, terminal local)

En otra ventana de terminal, **en tu laptop** (no en el clúster) ejecuta:

```bash
ssh -N -L 8787:node01:8794 TU_USUARIO@ken.lavis.unam.mx
```
> Hay algunos números en los que debes fijarte. A continuación qué es y de dónde lo obtienes, por si es necesario cambiarlos al ejecutar este comando

- `node01:8794` → nodo y puerto donde corre RStudio en el clúster (los que imprimió la Terminal 1).

> De nuevo, es importante no mover ni correr nada en esta segunda ventana, ambas las debes dejar abiertas para avanzar en los siguientes pasos.

Por último, en el navegador de tu computadora, escribir lo siguiente:
```text
http://localhost:8787
```

Con esto ya tienes RStudio corriendo en el clúster, con tus archivos y tu carpeta de trabajo.


### Paso 1. Ubicación de scripts y datos

> Debes analizar los siguientes documentos para poder tarbajar con estos datos:
'/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/regarding_data_route/data_as_route.md'
'/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/regarding_user_permissions/permissions.md'

Todos los scripts están en:

```text
/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts
```

> Nuevamente, copia la carpeta a tu directorio y trabaja ahí (no sobre la de Diana, no se debe ejecutar nada en ese directorio):

```bash
cp -r /mnt/data/bioinfo3/dbarrientos/proyecto_integrador /mnt/data/bioinfo3/TU_USUARIO/
```
> Los scripts tienen rutas absolutas a `/mnt/data/bioinfo3/dbarrientos/...` (`setwd(...)`, `file <- ...`). Después de copiar, edita esas rutas para que apunten a tu directorio.

### Paso 2. Librerías de R necesarias
Hay que descargar y cargar las siguientes librerías para poder trabajar. Para loomR y hdf5r, se explica abajo cómo se instala:
`hdf5r`, `Matrix`, `Seurat`, `SeuratDisk`, `rhdf5`, `SingleCellExperiment`, `loomR` (GitHub, rama `develop`), `scater`, `patchwork`, `devtools`/`remotes`.

#LoomR y hdf5r desde GitHub:
```r
install.packages("devtools")
devtools::install_github("hhoeflin/hdf5r")
devtools::install_github("mojaveazure/loomR", ref = "develop")
```
> Instala una sola vez. 
`loom_to_seurat_v2.r` trae esas instalaciones dentro del script. No las vuelvas a correr cada vez.

### Paso 3. Qué correr (orden recomendado)

1. `subset.r` → es el que **funciona** y produce un objeto Seurat de 10,000 células.
2. `read_seurat.r` → para leer el objeto recortado que dejó Yalbi.
3. `loom_to_seurat.r` → intento completo; hay que lanzarlo **como job** (ver su sección), no desde RStudio interactivo.
4. `v2`, `v3`, `v4` → solo como referencia de lo que se probó y por qué falló.
5. RECUERDA ANALIZAR PRIMERO 'QC_seurat_filtrado.qmd' Y DESPUÉS QC_seurat_yalbi_filtrado.qmd

| Dato | Dónde vive | Cómo se obtuvo |
|---|---|---|
| Objeto Seurat (`.rds`, versión recortada) | `/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-B.rds` | Preparado por la profesora y compartido en la carpeta común del curso |
| Archivo original (`.loom`, atlas completo) | `/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-B.loom` | Archivo original del dataset asignado, copiado a nuestra carpeta de trabajo |
| Objeto Seurat construido por nosotros (subset 10 k células) | /mnt/data/bioinfo3/dbarrientos/scripts
| Generado con `subset.r` |


## Conceptos clave del formato Loom (necesarios para entender los scripts)

Un archivo **Loom** es un contenedor HDF5 con esta estructura. La obtuvimos con `loomR` y su `help(loomR.package)`:
| Grupo / dataset | Tipo | Qué contiene |
|---|---|---|
| `matrix` | Dataset (float) | La matriz de expresión, **927,681 × 38,606** |
| `row_attrs` | Grupo | **Atributos de los renglones, de los genes/features** |
| `col_attrs` | Grupo | **Atributos de las columnas, de las células** |
| `row_graphs` | Grupo | Grafos entre genes (vecindades, redes). Aquí está vacío/no lo usamos |
| `col_graphs` | Grupo | Grafos entre células (como grafo de vecinos más cercanos). No lo usamos |
| `.matrix_dimnames` | Grupo | Nombres de dimensiones de la matriz |


> Otras aclaraciones:
- **`row_attrs`** = una tabla con una entrada por gen:
  - `ID`
  - `Symbol` → símbolo del gen (tiene duplicados y Seurat no acepta nombres de fila duplicados).
  - `Type` → tipo de feature.
  - `rownames` → IDs Ensembl, únicos (ejemplo: `ENSG00000290825`). Son 38,606.
- **`col_attrs`** = una tabla con una entrada por célula:
  - `colnames` → los barcodes de las células. Son 927,681
  - `colnames_factor` → los mismos barcodes como factor.
- **`row_graphs` / `col_graphs`** = grafos (conexiones entre genes / entre células), por ejemplo, qué células son vecinas. No los necesitamos para crear el objeto Seurat.


Equivalencias de Loom con Seurat:
| Loom | Seurat |
|---|---|
| Genes (`row_attrs`) | Filas = features |
| Células (`col_attrs`) | Columnas = cells (barcodes) |
| `matrix` (en nuestro archivo viene células × features | `counts` debe ir features × células, hay que transponerla |

## Estructura

| Carpeta / archivo | Qué va aquí |
|---|---|
| `reporte.qmd` | El reporte ejecutable (pendiente de crear) |
| `R/` | Funciones que uses en más de un lugar |
| `scripts/` | Los scripts descritos abajo (`loom_to_seurat*.r`, `subset.r`, `read_seurat.r`) |
| `datos/` | Vacía a propósito: los datos no se versionan |
| `figuras/` | Salida del reporte. Tampoco se versiona |
| `output/` | Salidas del script completo (`tabla_genes.rds`, objeto Seurat, `sessionInfo`) |

---


## Descripción de los scripts

### 1. `loom_to_seurat.r` — intento de conversión completa

Qué hace: intenta pasar todo el atlas (927,681 células) de Loom a Seurat leyendo la matriz por bloques. Es nuestro intento manual, guiado por (líneas ~104–110, no todas se usan):
<https://github.com/satijalab/seurat/blob/master/vignettes/conversion_vignette.Rmd>

> Se debe correr como job, ya que el archivo es demasiado pesado.
Se lanzó un job de clúster, `loom_to_seurat.sh`, que usa `module load apptainer/1.5.3`


**Paso a paso:**

1. Conectar al loom con `Connect()` en modo lectura. Ruta usada: `.../proyecto_integrador/data/proyecto-equipo-B.loom`.

2. Leer nombres: `row_attrs/Symbol` (símbolos), `row_attrs/rownames` (Ensembl, únicos) y `col_attrs/colnames` (barcodes). Imprime cuántos duplicados hay: los símbolos tienen duplicados; los Ensembl no. Por eso se usa Ensembl como nombre de fila del objeto.

3. Leer la matriz por bloques de 1,000 células (`tam_bloque <- 1000`):
   - Detecta la orientación (`celulas_en_filas`): si las células van en filas, transpone el bloque.
   - Cada bloque se convierte a matriz dispersa (`Matrix::Matrix(m, sparse = TRUE)`), que ocupa mucho menos memoria.
   - Con el primer bloque estima los no ceros totales y aborta si exceden ~2,147 millones, porque una `dgCMatrix` de R no admite más. En ese caso el mensaje pide "trabajar con un subconjunto de células".
   - Imprime el avance cada 10 bloques.
   - Une todos los bloques con `do.call(cbind, lista)`.

4. Tabla de equivalencias Ensembl: símbolo (`symbol_unico = make.unique(sym)` para resolver duplicados). Cierra el loom.

5. Guarda `tabla_genes.rds` en `output/` antes de seguir, por si algo falla después.

6. Crea el objeto Seurat con `CreateSeuratObject(counts = matriz_cuentas, project = "proyecto-equipo-B")` y agrega los símbolos como metadatos de genes dentro de un `tryCatch` (si falla solo avisa; la tabla ya está guardada).

7. Revisión rápida: `print`, `dim`, `head(tabla_genes)` y las primeras 5×5 de la matriz de cuentas.
**8. Guarda `proyecto-equipo-B_seurat.rds` y `sessionInfo` (en `output/`).**

Salidas esperadas (en `.../proyecto_integrador/output/`): 
`tabla_genes.rds`, `proyecto-equipo-B_seurat.rds`, `sessionInfo_loom_to_seurat.txt`.

> Pendiente: revisar el outlog del job `loom_to_seurat.sh` para ver dónde se detuvo. 


### 2. `loom_to_seurat_v2.r` — exploración con `loomR`

Este archivo sirvvió para encontrar los features → nombres de los genes (ensembl)
por ejemplo:

| Sección | Qué hace | ¿Funcionó? |
|---|---|---|
| Instalación de paquetes | `devtools::install_github` de `hdf5r` y `loomR` (rama `develop`) | Se hace una sola vez |
| `loomR::connect(..., skip.validate = TRUE)` | Abre el loom saltando la validación estricta del formato | sí |
| `names(lfilerow.attrs)`/`names(lfilecol.attrs)` | Lista los atributos: `"ID" "Symbol" "Type" "rownames"` y `"colnames" "colnames_factor"` | Como dijimos, así encontramos los features y los barcodes |
| Extraer `features` y `cells`, `table(Type)` | Revisa los tipos de feature antes de crear el Seurat | sí |
| `count_matrix <- lfile$matrix[, ]` → `t()` → `Matrix(sparse = TRUE)` → `CreateSeuratObject` | Intento de convertir todo de golpe | No se pudo (la matriz es demasiado grande) |
| Bloque `hdf5r` (`H5Filenew`,`f[["rowattrs"]]names`, `f[["col_attrs"]]$names`) | Inspeccionar qué contiene cada grupo y cómo se llama | Sí funcionó |
| `as.Seurat(lfile, assay = NULL, cells = "col_names", features = "ID")` | Conversión directa con `loomR` | No funcionó |
| "Crear subset" (`loom[["matrix"]][, 1:1000]`) | Leer solo 1,000 células para no cargar el atlas | Parcial |
| "Checar memoria" (`options(future.globals.maxSize = 16 * 1024^3)`, `Seurat.object.assay.version = "v5"`) | Subir el límite de memoria para objetos globales | sí ayuda, pero no resolvió |
| `SeuratDisk::LoadLoom(file, assay="RNA", cells="colnames", features="Symbol", filter="none")` | Carga directa con SeuratDisk | sí funcionó, problema de GB de la matriz |
| `as.loom(pbmc, ...)` y `as.Seurat(pbmc.loom)` + `VlnPlot` | Código de la viñeta de Seurat (dataset `pbmc`), no de nuestros datos | No se pudo crear `pbmc.loom` |

**Problemas de memoria que encontramos:** solo podíamos usar 1 de 32 CPUs y no toda la memoria disponible, así que se corrió algo para aumentar la memoria asignada. Aun así, la conversión completa no cabía.


### 3. `loom_to_seurat_v3.r` — loom reducido a 1,000 células

Qué hace:

1. Detecta los cores disponibles (`parallel::detectCores()`) y deja 2 libres para el sistema (`options(mc.cores = n_cores - 2)`). 

2. Abre el loom original con `loomR::connect(..., skip.validate = TRUE)`.

3. Recorta a 1,000 células (`n_use <- 1000`, `cells_use <- 1:n_use`) leyendo solo ese bloque: `mat_sub <- lfile[["matrix"]][cells_use, ]` (loomR presenta la matriz como células × genes).

4. Extrae metadatos de genes (`ID`, `Symbol`, `Type`, `rownames`) y de células (`colnames` solo de las 1,000 elegidas).

5. Crea un loom nuevo y estándar con `loomR::create(filename = "proyecto-equipo-B_1k.loom", data = mat_sub, gene.attrs = ..., cell.attrs = ..., do.transpose = FALSE, overwrite = TRUE)`.

6. Cierra los archivos.

> Logramos entender los metadatos de genes y de células, y se creó el loom reducido. **Salida:** `proyecto-equipo-B_1k.loom` (332 KB).

**Qué falló:** al agregar los metadatos de **células** aparece:

```text
Error in self$add.attribute(attributes = attributes, MARGIN = 2, overwrite = overwrite) :
  At least one dimmension for each cell attribute must be 38606
```

Lectura del error (hipótesis): `loomR` esperaba que los atributos de célula tuvieran longitud 38,606 (el número de genes), no 1,000. Es decir, interpretó la matriz con la orientación invertida: tomó las 38,606 columnas como si fueran células. Es un problema de orientación (`do.transpose = FALSE`), no de los metadatos en sí. Para probar: pasar `do.transpose = TRUE`, o transponer `mat_sub` antes de `create()`.


### 4. `loom_to_seurat_v4.r` — mismo intento, con diagnóstico

Nuevamente falló las dimensiones, sí se pudo crear el loom reducido, los metadatos se cargaron pero viene error:
Error in self$add.attribute(attributes = attributes, MARGIN = 2, overwrite = overwrite) : 
  At least one dimmension for each cell attribute must be 38606

**Qué hace:**
1. Configuración: detecta los cores (detectCores()) y deja 2 libres con options(mc.cores = n - 2). Hace setwd() a scripts/ y usa una ruta relativa para el loom ("proyecto-equipo-B.loom"), así que depende de que el setwd funcione y de que el .loom esté en esa carpeta.
2. Abre el loom original con loomR::connect(..., skip.validate = TRUE).
3. Elige el subset. Cuenta las células totales (length(lfile$col.attrs$colnames[]), 927,681) y usa n_use <- min(1000, n_total). Es una protección por si se pide más de lo que hay. Las células son las primeras 1,000 (seq_len).
4. Lee la matriz: mat_sub <- lfile[["matrix"]][cells_use, ]. Imprime dim(mat_sub) y el tamaño en RAM. Resultado: 1000 × 38606, o sea células × genes.
5. Arma metadatos mínimos: solo Symbol para genes y colnames (barcodes de las 1,000 células) para células.
6. Crea el loom reducido con loomR::create(..., do.transpose = FALSE, overwrite = TRUE, display.progress = TRUE). El nombre sale de n_use: proyecto-equipo-B_1000cells.loom.
7. Cierra y limpia: close_all(), rm() y gc().

Lo que cambió respecto a V3:
- En V3 la salida es: ‘_1k.loom’ y acá es: ‘_1000cells.loom’
- Los atributos ahora son solo “Symbol”
- Control de n_use
- Limpieza de objetos: rm() + gc()
- Ruta del loom es relativa, tras setwd.** A diferencia de V3 donde la ruta es absoluta**

print(dim(mat_sub)) muestra en la consola cuántas filas y cuántas columnas tiene la matriz mat_sub (la matriz que se extrajo del loom de las primeras 1,000 células)

Da como resultado:
Dimensión extraída:
	⁠print(dim(mat_sub)):
1000 células x 38606 genes


### 5. `subset.r` — el que funciona

Qué hace: crea un objeto Seurat con las primeras 10,000 células, sin necesidad de cargar el atlas completo. Usa `hdf5r` directamente (no `loomR`).

Paso a paso:

1. Configuración: `file` (ruta al `.loom`) y `n_use <- 10000L`.

2. Abrir el loom con `hdf5r::H5File$new(file, mode = "r")`.

3. Cómo se creó `d`:
   ```r
   mat_h5 <- h5[["matrix"]]   # dataset de expresión
   d <- mat_h5$dims           # dimensiones de la matriz
   ```
   `d` guarda las **dimensiones del dataset de expresión** (misma dimensión que la matriz original: 927,681 × 38,606). Se usa `$dims` y no `dim()` porque en `hdf5r` `dim()` no funciona sobre datasets HDF5.

4. Nombres: `genes` . Imprime cuántos hay (features y células).

5. Selección de células: `n_use. Las primeras 10,000, en el orden del archivo (no aleatorias).

6. Leer solo el subset, decidiendo por la orientación con `d`:
   - Si `d = features × células` → `mat_h5[, cells_use]`.
   - Si `d = células × features` → `mat_h5[cells_use, ]` y luego `t()`. Nuestro caso: células × features.)
   - Si no coincide con ninguna → `stop()` con mensaje de error.

7. Cerrar el HDF5 (`h5$close_all()`).

8. Convertir a dispersa (`Matrix::Matrix(counts, sparse = TRUE)`).

9. Asignar nombres con `make.unique()` (los símbolos duplicados quedan, por ejemplo, como `GEN`, `GEN.1`).

10. Resumen: `Matriz final: 38606 features x 10000 células`.

11. **Crear el Seurat:**
    ```r
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = counts, project = "proyecto-equipo-B",
      min.cells = 0, min.features = 0
    )
    ```
12. Guardar como **`proyecto-equipo-B_10000cells_seurat.rds`**.
13. QC mínima: un `VlnPlot` de `nFeature_RNA` y `nCount_RNA`, sin filtrar, solo para comprobar que el objeto sirve.

> El script usa una ruta relativa (`saveRDS(seurat_obj, output_file)` sin `setwd`), así que se guarda en el **directorio de trabajo de la sesión de R**. Por el listado de archivos, quedó en la carpeta `scripts/`.



### 6. `read_seurat.r` — leer el objeto recortado de la profesora

Qué hace: lee el `.rds` que preparó la profesora (la versión recortada del archivo original) y lo convierte a un objeto Seurat:

```r
data <- readRDS(".../scripts/proyecto-equipo-B.rds")
library(Seurat)
data_seurat <- as.Seurat(data, counts = "counts", data = NULL)
data_seurat
```

- counts toma la matriz de conteos del slot `counts`.
- `data = NULL` → no hay matriz normalizada; se normalizará después en Seurat.

Este script utiliza el trabajo de la profesora, mientras que  `loom_to_seurat.r` es nuestro intento manual con todo el dataset.

### 7. QC_seurat:
Dentro de estos documentos se hace el análisis final.
> ES IMPORTANTE QUE ESTE SEA EL PRIMERO QUE LEAS:
'/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/QC_seurat_filtrado.qmd'

> lUEGO SIGUE ESTE:
‘/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/QC_seurat_yalbi_filtrado.qmd’

###Archivos de salida en nuestro directorio:
| Archivo | Tamaño / fecha | Qué es | Generado por |
|---|---|---|---|
| `proyecto-equipo-B.loom` | 331 MB, 1 oct | Original: atlas completo (927,681 × 38,606) | Copiado de `compartido/10x/datasets/proyecto-2026/` |
| `proyecto-equipo-B.rds` | 17 MB, 2 oct | Objeto recortado de la profesora | Se lee con `read_seurat.r` |
| `proyecto-equipo-B_10000cells_seurat.rds` | 597 KB, 5 oct | Seurat con 10,000 células  | `subset.r` |
| `proyecto-equipo-B_1k.loom` | 332 KB, 5 oct 19:23 | Loom de 1,000 células (con el error de dimensiones) | `loom_to_seurat_v3.r` |
| `proyecto-equipo-B_1000cells.loom` | 332 KB, 5 oct 19:26 | Loom de 1,000 células (con diagnóstico) | `loom_to_seurat_v4.r` (por el nombre y la hora) |
| `proyecto-equipo-B_10k.loom` | 5 KB, 1 oct | Intento temprano de loom de 10 k células; casi vacío |  Ningún script que tengamos lo genera |

##Exploración inicial de los datos:

Debido a que tuvimos problemas con seurat disk, en donde no estaba instalado en R del cluster y en rstudio, donde estaba presente, había problemas de de espacio para convertirlo a seurat. Debido a esto, los profesores, en conjunto con lxs profesores, decidimos usar 3 opciones. 
1.- Hacer el analisis inicial con loom. Dado que se nos pide un archivo seurat para filtrar este archivo solo continuará hasta la exploración, a este lo llamaremos `loom` 
2.- Hacer una transformación seurat nosotros mismos usando rstudio en el cluster con un archivo cortado, a este lo llamaremos `cut` en esta seccion 
3.- Analizar el archivo transformado por los lxs doctores, el cuál fue igualmente cortado y transformado en rstudio, a este archivo lo llamaremos `yalbi` 

### Introducción: 
Para contestar estas preguntas, se utilizó el archivo loom proporcionado además de dos librerías de R previamente descargadas: `rhdf5` y `LoomExperiment`. El código en R para contestar las preguntas se encuentra en los scripts 

- **Para el archivo loom:** 
`/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/loom_analysis.R`,el cual se puede lanzar en SLURM con el job 
`/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/loom_analysis.sh`

- **Para el archivo cut:** 
`/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_analysis_cut.R`, el cual se puede lanzar en SLURRM con el job 
`/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_analysis_cut.sh`

- **Para el archivo Yalbi:** 
`/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_analysis_yalbi.R`, el cuál se puede lanzar en SLURM con el job 
`/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_analysis_yalbi.sh`

A continuación se explicará la lógica del código de cada sección y se explicarán los outputs para contestar las preguntas, los cuáles están disponibles en el la carpeta `/mnt/data/bioinfo3/dbarrientos/proyecto_integrador/out_logs/` de la forma:


- Para el archivo cut: seurat_analysis_cut_%j.log, seurat_analysis_cut_%j.err - Para el archvio Yalbi: seurat_analysis_yalbi_%j.log, seurat_analysis_yalbi_%j.err ### 1. ¿Son cuentas crudas?


En el documento previamente proporcionado, se especifica que si las cuentas poseen decimales, estas no son crudas dado que pudo haber una transformación o proceso de normalización. A parte de verificar si todos son decimales se verificará si todos los números son positivos, otro factor que puede indicar normalización o transformación (logaritmo, por ejemplo) 


Si revisamos el documento del de los logs para cada job, nos podemos dar cuenta que que en una sección dice: 


Para el archivo loom: 

```bash
Todos los conteos son no negativos y no enteros
```

Para el archivo seurat cut: 

```bash 
Todos los conteos son no negativos 
Todos los conteos son enteros 
```

Además de que en el summary nos arroja:

```bash 
Min. 1st Qu. Median Mean 3rd Qu. Max. 
0.0 0.0 1.0 26.6 4.0 11666.0 
N barcodes: 10000 | zeros: 3157 | min: 0 | max: 18483 
```

Para el archivo yalbi nos arroja: 

```bash 
Todos los conteos son no negativos 
Todos los conteos son enteros 
```

Y en el summary nos arroja: 

```bash 
Min. 1st Qu. Median Mean 3rd Qu. Max. 
0.00 0.00 1.00 29.99 4.00 13994.00 
N barcodes: 927681 | zeros: 290291 | min: 0 | max: 32191 
```

Todos estos son pruebas que para las 3 matrices, todos los conteos son crudos 

### 2. ¿Ya estan fitlradas las celulas? 
Para contestar estas preguntas, haremos gráficas de las cruvas de los UMI, en donde se ordenan los UMI y suman los UMIs totales de cada una de las ceuluas y se grafica una curva. Si es que observamos células con poca cantidad de transcritos totales, entonces podmemos concluir que no estan filtradas. 
Para el archivo loom, podemos hacerlo con emptydroplets directmanete, sin embargo para los archvios seurat, dado que estan cortados por la memoria, **emprtydroplets no permite hacer las graficas con tan pocos transcritos** pro lo tanto se desarrollo una funcion similar, 
que igualmente ordena y grafica los UMIs. Las graficas se pueden encotnrar en /mnt/data/bioinfo3/dbarrientos/proyecto_integrador/figures/ 

Si observamos en las 3, en todas hay células con pocos conteos, además, en el summary de cut y yalbi aparecia que había células con 0 transcritos, debido a esto podemos concluir que las células no están filtradas 

### 3. ¿Donde estan los metadados de lote, donante y condicion? 
Según los profesores durante el curso, no tiene sentido el buscar un lote o una condición debido a que no estamos haciendo analísis de expresión diferencial. Adicionalmente, al momento de consultar los metadatos en el archivo, se pudo apreciar que estos solo tenían referentes a los genes, no a las células. Para corroborar se incluyó en el código que dentro del archivo loom, se buscaran todas las palabras que hicieran alusión a estas caracterísiticas. 

Los resutlados para los 3 archivos fueron 


```bash 
character(0) 
```

Por lo cuál no se encotró ninguna información de lote o condicion, por lo cual se llegó a la conlcusión de que no había ninguna sección o bloque que se refirierar a metadatos, lote, donante o condición. Sin embargo pudimos llegar a la información del dataset original gracias al link: 
https://www.10xgenomics.com/datasets/peripheral-blood-mononuclear-cells-pbm-cs-from-a healthy-donor-manual-channel-5-3-1-standard-3-1-0 
Por lo cuál, en dado caso que se requiera información adicional, podemos consultarla aquí. 

### 4. ¿Qué versión de anotación de genes se usó? 
Consultando la información en la página web de los datos, podemos concluir que se utilizó la versión del genoma GRCh38-3.0.0 con la versión del pipeline de 3.1.0 y una química de Sinlge Cell 3` v3. Dado que los archivos seurat son submuestras del loom, los 3 tienen la misma anotación

## Preguntas iniciales 

Todas estas preugntas se pueden hacer con los out logs 
Para el archivo looom:
1. ¿Cuántos genes y cuántas células hay? genes: 38606 y células 927681 2. ¿Sobrevivieron los metadatos de lote? No. El archivo original no contiene metadatos de lote ni de donante. La verificación solo indica si la conversión preserva esas etiquetas, no dice nada biológico sobre lotes 
Para el archivo cut: 
1. ¿Cuántos genes y cuántas células hay? genes: 38606 y células 10000 
2. ¿Sobrevivieron los metadatos de lote? No. El archivo original no contiene metadatos de lote ni de donante. La verificación solo indica si la conversión preserva esas etiquetas, no dice nada biológico sobre lotes 
Para el archivo seurat Yalbi 
1. ¿Cuántos genes y cuántas células hay? genes: 38606 y células 927681 2. ¿Sobrevivieron los metadatos de lote? No. El archivo original no contiene metadatos de lote ni de donante. La verificación solo indica si la conversión preserva esas etiquetas, no dice nada biológico sobre lotes 
Aquí se hace la observacion que la mayoría son creos (290291)



## Bitácora de decisiones

Cada vez que tomes una decisión de análisis —un umbral, un método, un filtro—
anótala aquí con su razón. **La tabla de decisiones es el criterio A3 de la
rúbrica** (15 %), y llenarla al final, de memoria, se nota.

| Fecha | Decisión | Por qué | Qué se probó antes |
|---|---|---|---|
|21/Sep/2026 | Usar rstudio en lugar de correr como job |Las librerias necesarias para conversion de loom a seurat no estaban disponibles en R |Instalar las librerias, probablemente debido al espacio todos los intentos terminaron como "fatal error" |
|28/Sep/2026 | Tomar submuestras del conjunto de datos y convertirlas a seurat en lugar del dataset original ,una fué hecha por nosotros, otra por la Dra. Yalbi | El archivo loom completo no soportaba la conversion de loom a seurat en Rstudio, terminando como "fatal error" nuevamente | Hubieron muchos intentos de pasar de loom a seurat, inslusive algunos hechos por con ayuda de los profesores usando librerias para manejar la memoria |
|30/Sep/2026 | Utilizar un metodo prublicado por Lun et al. para estimar las gotas vacias en lugar de las metrias de knee e inflection | Debido a que teniamos pocas células, la librería de EmpryDrops fallaba al intentar estimar las celulas vacias | Se intento forzar e inclusive recrear las funciones de emprydrops usando su codigo fuente, sin embargo estos metodos no soportaban tan pocos valores |
|01/Oct/2026 | Mantener el nivel de significancia de 0.01 para el filtrado de celulas vacias usando emptydrops | Dado que estabamos usanod un metodo nuevo, ademas de que la distribucion de nuestros p-values lo sugeria, decidimos mantener este p-value | Graficar los p-values para comparar que se pierda, buscar articulos que isen el metodo de cuentiles a otro nivel menos estricto y perder menos datos|
|01/Oct/2026 | Aplicar filtros de QC para conservar células con más de 200 y menos de 2500 genes detectados (nFeature_RNA) y menos de 5 % de ARN mitocondrial (percent.mt) | Se establecieron estos umbrales para excluir células con pocas características detectadas como células de baja calidad, así como células con un número elevado de genes, que podrían corresponder a dobletes o multipletes, el límite de ARN mitocondrial permitió excluir células con una proporción elevada de transcritos mitocondriales, un posible indicador de daño celular o pérdida de integridad | Se evaluaron las métricas de control de calidad disponibles en el objeto Seurat, particularmente el número de genes detectados por célula y el porcentaje de ARN mitocondrial, para establecer criterios de filtrado antes de continuar con el análisis |
|03/Oct/2026 | No aplicar un filtro de QC basado en el porcentaje de ARN ribosomal (percent.ribo) | A diferencia del RNA mitocondrial, el porcentaje de ARN ribosomal no cuenta con un umbral de filtrado estandar, de hecho su proporcion ión puede variar según el tipo celular y las características de la muestra | Comparar el numero de celulas perdidas filtrando a traves de diferentes puntos de corte |
|05/Oct/2026 | Ajustar los parámetros de clustering y UMAP en el archivo cut para explorar la estructura del subset pequeño (k.param = 5, resolution = 1.5 y n.neighbors = 3) | Debido al tamaño reducido del subset (28 células), se modificaron los parámetros para explorar si era posible identificar y visualizar agrupaciones celulares, estos ajustes pueden producir agrupaciones artificiales, por lo que los clusters obtenidos no se consideran evidencia de poblaciones celulares reales ni resultados biológicos concluyentes | Se intentaron metrica para facilitar el clustering |
|05/Oct/2026 | Encontrar genes marcadores de cada clúster para el archivo cut, comparándolo contra todas las células restantes e informar únicamente los marcadores positivos | Identificar los genes que presentan mayor expresión en cada clúster para facilitar su caracterización y posterior identificación del tipo celular | Se utilizó una comparación de cada clúster contra el resto de las células |
|06/Oct/2026 | Seguir el flujo de análisis de scRNA-seq del tutorial de Seurat para el procesamiento, la normalización, la identificación de componentes principales, la construcción de grafos, el clustering y la visualización, realizar la anotación celular mediante marcadores génicos a partir de los clusters obtenidos | Se tomó como referencia el tutorial oficial de Seurat para el análisis de PBMCs, con el fin de mantener un flujo de trabajo reproducible y metodológicamente establecido | Se siguio el tutorial previamente establecido

> Referencias:
- Laboratorio de Biología Computacional (2026). Normalización, escalamiento y reducción de dimensionalidad. [Enlace](https://laboratorio-de-biologia-computacional.github.io/bioinfo-estadistica-3-lcg-2027-1-sitio/manual/assets/ceschin/clase-01/practica/Clase1_Normalizacion_y_Reduccion_de_Dimensionalidad.html)
- Satija Lab (apr, 28 - 2026). Seurat - Guided Clustering Tutorial. [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
