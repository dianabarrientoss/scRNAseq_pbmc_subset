
---
title: "COMANDO PARA LA RUTA DE NUESTROS DATOS"
output: html_notebook
author: "Barrientos-González Diana, Cervantes-Levario Román & Contreras-Rivera Jair"
date: "20-09-2026"
---

La ruta apunta a un directorio del clúster que solo tiene sentido dentro de ken, el cual tiene como objeto nuestros datos

```{bash}
(base) [dbarrientos@node01 data]$ ln -s /mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-B.loom data
(base) [dbarrientos@node01 data]$ ls -lh
total 0
lrwxrwxrwx 1 dbarrientos bioinfo3 79 Sep 20 22:28 data -> /mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-B.loom
```

Los archivos seurat se crearon a partir del archivo loom, esta las podemos encontrar dentro de  scripts. Es necesario copiarlos a su carpeta personal debido a que es muy costoso replicar la transformacion a .rds 

```{bash}
#Cut
cp /mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/proyecto-equipo-B_10000cells_seurat.rds /mnt/data/bioinfo3/usuario/proyecto_integrador/scripts

cp /mnt/data/bioinfo3/dbarrientos/proyecto_integrador/scripts/seurat_yalbi.rds /mnt/data/bioinfo3/usuario/proyecto_integrador/scripts
```
