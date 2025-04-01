# PIPELINE ANÁLISIS DE EXPRESIÓN DIFERENCIAL COMPLETO
En esta página se detalla un pipeline completo de análisis de datos de RNA-Seq (incluyendo descarga de muestras de GEO, pre-procesamiento de datos crudos, control de calidad, alineamiento, cuantificación y análisis de expresión diferencial) empleado para la realización del TFM titulado *"Estudio transcriptómico en epitelio nasal del efecto del sexo en el asma"*.
Los scrips necesarios se encuentran en la carpeta scripts de este repositorio. También es posible descargar el entorno de conda con todos los programas empleados y en el archivo README encontrarán el código para instalar y cargar todos los paquetes de R que se emplearán en este pipeline.
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
