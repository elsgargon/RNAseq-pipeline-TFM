# RNAseq-pipeline-TFM
En este repositorio se incluyen todos los script empleados para el desarrollo del TFM "Estudio transcriptómico" en el máster en Bioinformática de la Universidad VIU. 
Los scripts se emplean para hacer un análisis de expresión diferencial entre asmáticos y controles estratificado por sexo. El orden de ejecución de scripts viene determinado por su nombre "1_Nombre". 

## **1. Programas empleados en el preprocesado**:
  - **SRAToolkits**: para la descarga de archivos de Gene Expression Omnibus
  - **FastQC**: para el control de calidad de las lecturas.
  - **MultiQC**: para la visualización conjunta de los ficheros de control de calidad de todas las muestras.
  - **FastP**: para estimar el número de lecturas que presentan una base indeterminada (N).
  - **Fastq-Screen**: para determinar la existencia de posibles contaminaciones de otros organismos.
  - **Trimmomatic**: para realizar el filtrado de las lecturas por calidad.
  - **STAR**: para alinear las lecturas al genoma de referencia.
  - **Bamtools**: para realizar un control de calidad de las lecturas alineadas (ficheros .bam).
  - **RSEM**: para cuantificar la expresión de los genes.
  - **RNASeQC**: ?

Todos estos pogramas empleados se encuentran ya descargados en el entorno de conda RNAseq_pipeline.yml para mayor reproducibilidad. Descarga el archivo e introduce el siguiente comando en el terminal para crear el entorno de conda desde el archivo:
```markdown
conda env create -f RNAseq_pipeline.yml
conda activate RNAseq_pipeline


## **2. Paquetes de R empleados en el análisis de expresión diferencial** 

Deseq2, tximport, sva, readr, pheatmap, ggplot2, IHW, AnnotationDbi, org.Hs.eg.db, EnhancedVolcano, OUTRIDER, clusterProfiler, enrichplot y tidyverse.

El siguiente código comprueba si todos los paquetes necesarios para correr el script *5_Analsis_expresion_diferencial* en R están descargados. En caso de faltar alguno, se descargará automáticamente. Además, este código carga todos los paquetes necesarios en memoria, quedando todo listo para ejecutar el  script.
```r
###################################################################################
# Comprobación, descarga y carga en memoria automática de los paquetes necesarios #
###################################################################################

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


