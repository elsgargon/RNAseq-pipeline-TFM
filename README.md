# RNAseq-pipeline-TFM
En este repositorio se incluyen todos los script empleados para el desarrollo del TFM "Estudio transcriptómico" en el máster en Bioinformática de la Universidad VIU, edición 2024-2025. 
Los scripts se emplean para hacer un análisis de expresión diferencial entre asmáticos y controles estratificado por sexo. El orden de ejecución de scripts viene determinado en su nombre. En el apartado PIPELINE.md se encuentran detallados los pasos necesarios para ejecutar el pipeline al completo. 

A continuación encontrarás los pasos previos necesarios para la descarga de los programas y paquetes de R empleados en este pipeline.

## **1. Programas empleados en el preprocesado**:

<table align="center">
  <tr>
    <th>Herramienta</th>
    <th>Uso en el pipeline</th>
  </tr>
  <tr>
    <td>SRAToolkits</td>
    <td>Para la descarga de archivos de Gene Expression Omnibus</td>
  </tr>
  <tr>
    <td>FastQC</td>
    <td>Para el control de calidad de las lecturas.</td>
  </tr>
  <tr>
    <td>MultiQC</td>
    <td>Para la visualización conjunta de los ficheros de control de calidad de todas las muestras.</td>
  </tr>
  <tr>
    <td>FastP</td>
    <td>Para estimar el número de lecturas que presentan una base indeterminada (N).</td>
  </tr>
  <tr>
    <td>Fastq-Screen</td>
    <td>Para determinar la existencia de posibles contaminaciones de otros organismos.</td>
  </tr>
  <tr>
    <td>Trimmomatic</td>
    <td>Para realizar el filtrado de las lecturas por calidad.</td>
  </tr>
  <tr>
    <td>STAR</td>
    <td>Para alinear las lecturas al genoma de referencia.</td>
  </tr>
  <tr>
    <td>Bamtools</td>
    <td>Para realizar un control de calidad de las lecturas alineadas (ficheros .bam).</td>
  </tr>
  <tr>
    <td>RSEM</td>
    <td>Para cuantificar la expresión de los genes.</td>
  </tr>
  <tr>
    <td>RNASeQC</td>
    <td>?</td>
  </tr>
</table>

Todos estos pogramas empleados se encuentran ya descargados en el entorno de conda RNAseq_pipeline.yml para mayor reproducibilidad. Asegúrate de tener Conda instalado. Descarga el archivo e introduce el siguiente comando en el terminal para crear el entorno de conda desde el archivo:
```markdown
conda env create -f RNAseq_pipeline.yml
conda activate RNAseq_pipeline
```

## **2. Paquetes de R empleados en el análisis de expresión diferencial** 

Los paquetes empleados son:

<table align="center">
  <tr>
    <td>DESeq2</td>
    <td>AnnotationDbi</td>
    <td>tximport</td>
    <td>org.Hs.eg.db</td>
    <td>sva</td>
  </tr>
  <tr>
    <td>EnhancedVolcano</td>
    <td>readr</td>
    <td>OUTRIDER</td>
    <td>pheatmap</td>
    <td>clusterProfiler</td>
  </tr>
  <tr>
    <td>ggplot2</td>
    <td>enrichplot</td>
    <td>IHW</td>
    <td>tidyverse</td>
  </tr>
</table>

El siguiente código comprueba si todos los paquetes necesarios para correr el script *5_Analisis_expresion_diferencial* en R están descargados. En caso de faltar alguno, se descargará automáticamente. Además, este código carga todos los paquetes necesarios en memoria, permitiendo ejecutar el script sin hacer ningún paso adicional.
```r
#Lista de paquetes necesarios
required_packages <- c("DESeq2", "tximport", "sva", "readr", "pheatmap", "ggplot2", 
                       "IHW", "AnnotationDbi", "org.Hs.eg.db", "EnhancedVolcano", 
                       "OUTRIDER", "clusterProfiler", "enrichplot", "tidyverse")

#Función para instalar paquetes faltantes
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Instalando paquetes faltantes: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  } else {
    message("Todos los paquetes están instalados.")
  }
}

#Instalar paquetes de Bioconductor si faltan
install_bioconductor_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Instalando paquetes de Bioconductor faltantes: ", paste(missing_packages, collapse = ", "))
    BiocManager::install(missing_packages, update = FALSE, ask = FALSE)
  }
}

#Paquetes de CRAN
cran_packages <- c("readr", "pheatmap", "ggplot2", "IHW", "tidyverse", "enrichplot")

#Paquetes de Bioconductor
bioconductor_packages <- c("DESeq2", "tximport", "sva", "AnnotationDbi", "org.Hs.eg.db", 
                           "EnhancedVolcano", "OUTRIDER", "clusterProfiler")

#Instalar paquetes de CRAN
install_if_missing(cran_packages)

#Instalar paquetes de Bioconductor
install_bioconductor_if_missing(bioconductor_packages)

#Cargar paquetes en memoria
lapply(required_packages, library, character.only = TRUE)


