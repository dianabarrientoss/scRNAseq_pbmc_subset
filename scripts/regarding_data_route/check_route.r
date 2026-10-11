



ruta_datos <- "/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-B.loom"

if (!nzchar(ruta_datos)) {
  cat("Sin ruta: la practica va a usar el objeto SINTETICO de respaldo.\n",
      "Sirve para aprender el manejo de los objetos, pero las cifras no son\n",
      "las del experimento. Pega la ruta del chat para trabajar con los datos.\n")
} else if (file.exists(ruta_datos)) {
  cat("Listo: el objeto existe y se puede leer.\n")
} else {
  cat("La ruta esta escrita pero NO existe. Revisa que este completa y\n",
      "que estes dentro de ken. Avisa en el chat antes de seguir.\n")
}
