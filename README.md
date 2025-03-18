# RNAseq-pipeline-TFM
En este repositorio se incluyen todos los script empleados para el desarrollo del TFM "Estudio transcriptómico" en el máster en Bioinformática de la Universidad VIU. 
Los scripts se emplean para hacer un análisis de expresión diferencial entre asmáticos y controles estratificado por sexo. El orden de ejecución de scripts viene determinado por su nombre "1_Nombre". 

Se emplean los siguientes **programas**:
  - **SRAToolkits**: para la descarga de archivos de Gene Expression Omnibus
  - **FastQC**: para el control de calidad de las lecturas.
  - **MultiQC**: para la visualización conjunta de los ficheros de control de calidad de todas las muestras.
  - **FastP**: para estimar el número de lecturas que presentan una base indeterminada (N).
  - **Trimmomatic**: para realizar el filtrado de las lecturas por calidad.
  - **STAR**: para alinear las lecturas al genoma de referencia.
  - **Bamtools**: para realizar un control de calidad de las lecturas alineadas (ficheros .bam).
  - **RSEM**: para cuantificar la expresión de los genes.
  - **RNASeQC**: ?
  - **Deseq2**: para el análisis de expresión diferencial en R. 

**Paquetes de R** empleados en el análisis de expresión diferencial: 
- tximport
- sva
- readr
- pheatmap
- ggplot2
- IHW
- AnnotationDbi
- org.Hs.eg.db
- EnhancedVolcano
- OUTRIDER
- clusterProfiler
- enrichplot
- tidyverse

