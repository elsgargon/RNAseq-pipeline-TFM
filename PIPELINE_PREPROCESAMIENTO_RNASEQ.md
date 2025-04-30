# PIPELINE PREPROCESAMIENTO DE DATOS DE RNASEQ
En esta página se detalla un pipeline completo (incluyendo descarga de muestras de GEO, pre-procesamiento de datos crudos, control de calidad, alineamiento y cuantificación) empleado para la realización del TFM titulado *"Estudio transcriptómico en epitelio nasal del efecto del sexo en el asma"*.
Los scrips necesarios se encuentran en la carpeta scripts de este repositorio. También es posible descargar el entorno de conda con todos los programas empleados y en el archivo README encontrarán el código para instalar y cargar todos los paquetes de R que se emplearán.
## Estructura final de directorios
```plaintext
.
├── data
│   ├── raw
│   ├── processed
│   │   ├── Control_calidad
│   │   │   ├── multiqc
│   │   │   ├── contaminaciones
│   │   │   │   └── multiqc
│   │   │   ├── fastp_N_reads
│   │   │   │   └── multiqc
│   │   ├── Alineamiento
│   │   │   └── Control_calidad_bam
│   │   └── cuantificacion
│   │   └── discordancias_sexo
│   └── reference
└── scripts
│   ├── 1_descarga_y_conversión.sh
│   ├── 2_control_calidad.sh
│   ├── 3_alineamiento_control_calidad_bam.sh
│   └── 4_cuantificacion_rsem.sh
└── entornos_conda
    ├── RNAseq_pipeline.yml
    ├── FastQScreen.yml
    └── RSEM.yml
```

## 1. Descarga de muestras desde GEO y conversión a fastq
Es necesario **descargar previamente el "Accesion List" de GEO** y almacenarla en el **directorio data/raw** bajo el nombre de **muestras.txt**. Los archivos resultantes se almacenarán en esa misma carpeta en formato .fastq.gz para mejor optimización del espacio de almacenamiento. Para descargar y convertir las muestras de forma automática, ejecutar:
```markdown
./scripts/1_descarga_y_conversión.sh
```
## 2. Control de calidad
Analizamos la calidad de las lecturas crudas descargadas empleando diferentes programas:
- **FastQC**: control de calidad general de lecturas de NGS.
- **Fastqscreen**: comprobar si existen contaminaciones en nuestra muestra usando los genomas predeterminados.
- **fastp**: en nuestro caso existía un alto número de bases indeterminadas en algunas posiciones (N) por lo que empleamos este programa para contabilizar el número total de N de las lecturas.
- **MultiQC**: visualización conjunta de los resultados de todas las muestras.

Estos informes deben revisarse para evaluar si es necesario realizar algun filtrado o recorte de las lecturas para mejorar su calidad. En nuestro caso no fue necesario.
Para realizar este control de calidad de forma automática ejecutar:
```markdown
./scripts/2_control_calidad.sh
```
## 3. Alineamiento al genoma de referencia y control de calidad del alineamiento
Existen diferentes herramientas para realizar el alineamiento de las lecturas. En este caso empleamos el programa STAR junto con el genoma de referencia **GRCh38.p14.genome.fa** y el archivo de anotación **gencode.v47.annotation.gtf** para el alineamiento. Es necesario descargar estos dos archivos previamente y almacenarlos en la carpeta data/reference. Bamtools se emplea para el control de calidad de los archivos .bam generados.
Para realizar el alineamiento y el control de calidad de los archivos .bam, ejeutar:
```markdown
./scripts/3_alineamiento_control_calidad_bam.sh
```

En este caso realizamos un paso adicional no incluido en el script automátizado en el que comprobamos si las muestras que presentan un alto número de indeterminaciones presentan asociación significativa con la calidad del mapeo (podría indicar la necesidad de filtrar o recortar estas bases si se relacionara con una calidad de mapeo baja): 
```markdown
#####################################################################################
# Confirmar que las indeterminaciones (N) no están afectando al porcentaje de mapeo #
#####################################################################################

# Extraer % de mapeo por muestra
folder_path="data/processed/Alineamiento/Control_calidad_bam" # Directorio donde se encuentran los archivos de salida de Bamtools
output_file="data/processed/Alineamiento/Control_calidad_bam/porcentaje_mapeo.txt" # Fichero donde se almacenarán los porcentajes de mapeo de cada muestra

# Encabezado del archivo de salida
echo -e "SRR_ID\tMapped_Percentage" > "$output_file"
# Iterar sobre los archivos .txt en la carpeta
for file in "$folder_path"/*.txt; do
    srr_id=$(basename "$file" | sed -E 's/^(SRR[0-9]+)_.*$/\1/')
    mapped_percentage=$(grep "Mapped reads" "$file" | sed -E 's/.*\(([0-9]+\.[0-9]+)%\).*/\1/')
    echo -e "$srr_id\t$mapped_percentage" >> "$output_file"
done
```
```r
########
# EN R #
########
# Añadir una columna manualmente con excel llamada n_content_fail_warn que indique que muestras presentan un "Warn" o "Fail" en el contenido de N del reporte de FastQC. Poner 1 si lo tienen y 0 si no.

# Cargar datos
data <- read.table("data/processed/Alineamiento/Control_calidad_bam/porcentaje_mapeo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Hacer regresión GLM
modelo <- glm(Mapped_Percentage ~ n_content_fail_warn, data = data)
summary(modelo)
```
## 4. Cuantificación de la expresión génica
Una vez alineadas las lecturas, podemos realizar la cuantificación de la expresión génica, tanto a nivel de gen como a nivel de exones, con el programa RSEM. Los resultados se almacenarán en una nueva carpeta llamada "cuantificacion" Para realizar la cuantificación de forma automática ejecutar:
```markdown
./scripts/4_cuantificacion_rsem.sh
```
## 5. Detección de discordancias de sexo
Dado que vamos a realizar un análisis estratificado por sexo, comprobaremos que el sexo informado por el paciente corresponde al sexo inferido mediante los datos de expresión. Para ello emplearemos la expresión del gen XIST (cromosoma X) y de los genes del cromosoma Y EIF1AY, KDM5D, UTY, DDX3Y y RPS4Y1. 

Clasificaremos como mujeres aquellas muestras que presenten una expresión del gen XIST al menos dos veces mayor a la suma de la expresión de los genes del cromosoma Y. Clasificaremos como hombres aquellas muestras cuya suma de la expresión de los genes del cromosoma Y sean al menos dos veces mayor a la expresión de XIST. Si no se da ninguno de estos dos casos, la muestra se clasifica como indeterminada. Los individuos que presenten un discordancia entre el sexo informado y el calculado deben ser eliminados. 

Necesitaremos crear un archivo llamado sexo_informado.txt que incluya el nombre de las muestras y el sexo informado, en el mismo formato que el archivo "sexo_informado.txt" que se encuentra en la carpeta data/raw de este repositorio.

Una vez hecho esto, emplearemos el siguiente código para inferir el sexo de los individuos e identificar discordancias con el sexo informado:
```r
# 1. Ruta de la carpeta donde están los archivos de RSEM
folder_path <- "data/processed/cuantificación"

# 2. Obtener la lista de archivos .genes.results en la carpeta
file_list <- list.files(folder_path, pattern = "*.genes.results", full.names = TRUE)

# 3. Inicializar una lista para guardar los resultados
results <- list()

# 4. Leer el archivo de sexo informado
sexo_informado <- read.table("data/raw/sexo_informado.txt", sep = "\t", header = TRUE)

# 5. Leer el archivo de mapeo de Gene ID a Gene Name
gene_mapping <- read_csv("data/reference/ensemble_id_to_gene_name.txt")

# 6. Procesar cada archivo
for (file_path in file_list) {
  
  # 7. Cargar los datos de expresión génica
  gene_expression <- read_delim(file_path, delim = "\t")
  
  # 8. Realizar la correspondencia de IDs de genes con los nombres de genes
  gene_expression <- merge(gene_expression, gene_mapping, by.x = "gene_id", by.y = "Gene_stable_ID_version")
  
  # 9. Seleccionar solo los genes de interés
  genes_interes <- c("XIST", "EIF1AY", "KDM5D", "UTY", "DDX3Y", "RPS4Y1")
  expr_genes_interes <- gene_expression %>%
    filter(Gene_name %in% genes_interes)
  
  # 10. Convertir los datos de expresión a un formato que facilite el análisis
  expr_genes_interes <- expr_genes_interes %>%
    select(Gene_name, FPKM)
  
  # 11. Extraer el ID de la muestra desde el nombre del archivo (por ejemplo, SRR11951898)
  sample_id <- str_extract(file_path, "SRR\\d+")
  
  # 12. Filtrar para obtener el sexo informado para la muestra correspondiente
  sexo_informado_muestra <- sexo_informado %>%
    filter(Muestras == sample_id) %>%
    pull(Gender)
  
  # 13. Pivotar los datos para obtener un formato adecuado
  expr_genes_interes <- expr_genes_interes %>% pivot_wider(names_from = Gene_name, values_from = FPKM)
  
  # 14. Calcular el sexo basado en la expresión génica
# Calcular la suma de los genes del cromosoma Y
suma_genes_Y <- rowSums(expr_genes_interes[, c("EIF1AY", "KDM5D", "UTY", "DDX3Y", "RPS4Y1")], na.rm = TRUE)

  sexo_calculado <- ifelse(
  expr_genes_interes$XIST >= 2 * suma_genes_Y, "Female",
  ifelse(suma_genes_Y >= 2 * expr_genes_interes$XIST, "Male", "Indeterminado")
  )
  
  # 15. Comparar el sexo calculado con el sexo informado y mostrar el resultado
  discordancia <- ifelse(sexo_calculado != sexo_informado_muestra, "Sí", "No")
  
  # 16. Almacenar los resultados en la lista
  results[[sample_id]] <- list(
    expr_data = expr_genes_interes,
    sexo_informado = sexo_informado_muestra,
    sexo_calculado = sexo_calculado,
    discordancia = discordancia
  )
}

# 17. Crear una lista para almacenar los resultados combinados
combined_results <- list()

# 18. Iterar sobre cada elemento en 'results' y extraer la información necesaria
for (sample_id in names(results)) {
  
  # Obtener los datos de expresión genética para cada muestra
  expr_data <- results[[sample_id]]$expr_data
  
  # Crear una fila con los valores de expresión y los metadatos (sexo calculado y discordancia)
  expr_data_with_metadata <- expr_data %>%
    mutate(Sample = sample_id, 
           Sexo_Calculado = results[[sample_id]]$sexo_calculado,
	   Sexo_Informado= results[[sample_id]]$sexo_informado,
           Discordancia = ifelse(length(results[[sample_id]]$discordancia) > 0, results[[sample_id]]$discordancia, NA))
  
  # Añadir la fila a la lista de resultados combinados
  combined_results[[sample_id]] <- expr_data_with_metadata
}

# 19. Convertir la lista de resultados en un data frame
final_results <- bind_rows(combined_results)

# 20. Reorganizar las columnas para que 'Sample' sea la primera columna
final_results <- final_results %>%
  select(Sample, everything())  # 'everything()' mantiene el resto de las columnas tal como están

# 22. Guardar los resultados en un archivo CSV
write.table(final_results, "data/processed/discordancias_sexo/sexo_inferido_new_rnaseq.csv", quote=F, row.names=F)

# 23. Mostrar solo las muestras con discordancia
muestras_discordantes <- final_results %>%
  filter(Sexo_Calculado != Sexo_Informado)

if (nrow(muestras_discordantes) > 0) {
  cat("Se encontraron muestras con discordancia entre el sexo informado y el calculado:\n")
  print(muestras_discordantes$Sample)
} else {
  cat("No se encontraron discordancias entre el sexo informado y el calculado.\n")
}

# 24. Representación gráfica del resultado
# Cambiar el formato de 'final_results' a largo para la gráfica de barras apiladas
data_long <- final_results %>%
  select(Sample, XIST, EIF1AY, KDM5D, UTY, DDX3Y, RPS4Y1, Sexo_Calculado, Sexo_Informado) %>%
  pivot_longer(cols = XIST:RPS4Y1, names_to = "Gen", values_to = "Expresion_TPM")

# Calcular la posición en y máxima para las etiquetas de cada muestra
max_y_position <- data_long %>%
  group_by(Sample) %>%
  summarize(Max_Expresion = sum(Expresion_TPM, na.rm = TRUE))

# Crear la gráfica y guardarla como PNG
output_path <- "data/processed/discordancias_sexo/sexo_muestras_calculado_grafica.png"

# Crear una columna con el sexo calculado y un asterisco si hay discordancia
final_results <- final_results %>%
  mutate(Sexo_Label = ifelse(Sexo_Calculado != Sexo_Informado,
                             paste0(Sexo_Calculado, "*"),
                             Sexo_Calculado))

annot_data <- max_y_position %>%
  mutate(Label_Pos = Max_Expresion * 1.1) %>%  # 10% por encima del máximo
  left_join(final_results %>% select(Sample, Sexo_Label), by = "Sample")

# Crear el gráfico
ggplot(data_long, aes(x = Sample, y = Expresion_TPM, fill = Gen)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("XIST" = "skyblue", "EIF1AY" = "orange", 
                               "KDM5D" = "purple", "UTY" = "green", 
                               "DDX3Y" = "pink", "RPS4Y1" = "yellow")) +
  labs(x = "Muestras", y = "Expresión Génica (TPM)", title = "Expresión de genes de cromosomas sexuales por muestra") +
  theme_classic() +
  # Etiquetas de sexo calculado (con * si hay discordancia)
  geom_text(data = annot_data,
            aes(x = Sample, y = Label_Pos, label = Sexo_Label),
            inherit.aes = FALSE, color = "black", size = 3, angle = 90) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Guardar la gráfica como archivo PNG
ggsave(output_path, width = 10, height = 6)

```
## 6. Análisis de expresión diferencial estratificado por sexo con Deseq2
Realizaremos en análisis estadístico para determinar si existen genes diferencialmente expresados entre hombre y mujeres con y sin asma. Para ello emplearemos el paquete de R Deseq2 basándonos el pipeline de bioconductor (https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow), aunque se incluyen ciertos pasos adicionales.

```r

```
