# PIPELINE ANÁLISIS DE EXPRESIÓN DIFERENCIAL COMPLETO
En esta página se detalla un pipeline completo de análisis de datos de RNA-Seq (incluyendo descarga de muestras de GEO, pre-procesamiento de datos crudos, control de calidad, alineamiento, cuantificación y análisis de expresión diferencial) empleado para la realización del TFM titulado *"Estudio transcriptómico en epitelio nasal del efecto del sexo en el asma"*.
Los scrips necesarios se encuentran en la carpeta scripts de este repositorio. También es posible descargar el entorno de conda con todos los programas empleados y en el archivo README encontrarán el código para instalar y cargar todos los paquetes de R que se emplearán en este pipeline.
## 1. Descarga de muestras desde GEO y conversión a fastq
Es necesario **descargar previamente el "Accesion List" de GEO** y almacenarla en el **directorio data/raw** bajo el nombre de **muestras.txt**. Los archivos resultantes se almacenarán en esa misma carpeta en formato .fastq.gz para mejor optimización del espacio de almacenamiento. Para descargar y convertir las muestras de forma automática, ejecutar:
```markdown
./scripts/1_descarga_y_conversión
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
./scripts/2_control_calidad
```
## 3. Alineamiento al genoma de referencia
Existen diferentes herramientas para realizar el alineamiento de las lecturas. En este caso empleamos el programa STAR junto con el genoma de referencia GRCh38.p14 y el archivo de anotación genómica
