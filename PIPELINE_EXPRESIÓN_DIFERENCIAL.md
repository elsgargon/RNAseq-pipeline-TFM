# PIPELINE ANÁLISIS DE EXPRESIÓN DIFERENCIAL ESTRATIFICADO POR SEXOS
En esta página se detalla un pipeline completo de análisis de expresión diferencial estratificado por sexos usando Deseq2. Se ha empleado para la realización del TFM titulado *"Estudio transcriptómico en epitelio nasal del efecto del sexo en el asma"*. Se detalla tanto el script de análisis de mujeres como de hombres por separado.
Los scrips necesarios se encuentran en la carpeta scripts de este repositorio. También es posible descargar el entorno de conda con todos los programas empleados y en el archivo README encontrarán el código para instalar y cargar todos los paquetes de R que se emplearán.

## Análisis de expresión diferencial con Deseq2
Realizaremos en análisis estadístico para determinar si existen genes diferencialmente expresados entre hombre y mujeres con y sin asma. Para ello emplearemos el paquete de R Deseq2 basándonos el pipeline de bioconductor (https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow), aunque se incluyen ciertos pasos adicionales.

## Mujeres
```r
#####################
# 1. Importar datos #
#####################

# Crear archivo con rutas de archivos rsem
while read identifier; do
  find data/processed/cuantificacion/ -type f -name "${identifier}_.isoforms.results" -exec cp {} data/processed/cuantificacion/mujer/ \;
done < women_id_newrnaseq.txt

# Establecer directorio
setwd("data/processed/cuantificacion/mujer/")

# Cargar archivos y establecer rutas
dir <- "data/processed/cuantificacion/mujer/"
samples <- read.table("data/raw/metadata_new_rnaseq_sinbajacalidad.txt", header = TRUE)

# Cambiar formato de datos a binario y filtrar mujeres
samples$Gender <- ifelse(samples$Gender == "Female", 1, ifelse(samples$Gender == "Male", 0, samples$Gender))
samples<-samples[samples$Gender==1, ]

# Convertir variables a factores
samples$Status <- factor(samples$Status)
samples$Gender <- factor(samples$Gender)
samples$Ethnicity <- factor(samples$Ethnicity)
samples$current_smoker <- factor(samples$current_smoker)
samples$Age <- factor(samples$Age)
samples$Smoke_Ever <- factor(samples$Smoke_Ever)
samples$smoke_pack_years_1 <- factor(samples$smoke_pack_years_1)

# Poner nombres de muestra 
rownames(samples) <- samples$V1.y
files <- paste0(dir, "/", samples$V1.y, "_.genes.results")
names(files) <- samples$V1.y

# Realizar test de correlación con variable Status para seleccionar covariables (aquellas significativas deben incluirse en el modelo)
cor.test(as.numeric(samples$Status), as.numeric(samples$Age))   
cor.test(as.numeric(samples$Status), as.numeric(samples$Ethnicity)) 
cor.test(as.numeric(samples$Status), as.numeric(samples$BMI)) 
cor.test(as.numeric(samples$Status), as.numeric(samples$current_smoker)) 
cor.test(as.numeric(samples$Status), as.numeric(samples$Smoke_Ever)) 
cor.test(as.numeric(samples$Status), as.numeric(samples$smoke_pack_years_1)) 
cor.test(as.numeric(samples$Status), as.numeric(samples$total_IgE))

# Create tx2gene mapping table a partir de los datos de RSEM
tx2gene <- data.frame(
  TXNAME = rownames(read.delim(files[1], nrows = 1, header = TRUE)),
  GENEID = rownames(read.delim(files[1], nrows = 1, header = TRUE))
)

# Importar datos con tximport
txi <- tximport(files, type = "rsem", tx2gene = tx2gene)

# Identificar genes con longitud 0 en todas las muestras
genes_cero <- rowSums(txi$length == 0) == ncol(txi$length)
sum(genes_cero)

# Recomendado por el creador, sustituir 0 por 1 (https://www.biostars.org/p/293639/) para permitir la carga de los datos
txi$length[txi$length == 0] <- 1

# Crear DESeqDataSet con DESeq2
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Status + Age + Ethnicity + BMI) 

# Calcular variantes surrogadas (SVs) con SVAseq 
# 1. Eliminar genes con bajo conteo (se recomienda especificar el tamaño del grupo más pequeño como el número mínimo de muestras).
smallestGroupSize <- 9
keep <- rowSums(counts(ddsTxi) >= 10) >= smallestGroupSize
ddsTxi <- ddsTxi[keep,]

# 2. Normalización y transformación de los datos para calcular SVs
ddsTxi <- estimateSizeFactors(ddsTxi)
norm_counts <- counts(ddsTxi, normalized=TRUE)

# 3. Calcular variantes surrogadas
mod <- model.matrix(~ Status + Age + Ethnicity, samples)    # Modelo completo
mod0 <- model.matrix(~ Age + Ethnicity, samples)    # Modelo nulo
svseq_1 <- svaseq(norm_counts, mod, mod0)
                    
# Añadir SVs al diseño experimental de DESeq2 y volver a crear DESeqDataSet con DESeq2
samples <- cbind(samples, svseq_1$sv)
names(samples)[(ncol(samples)-1):ncol(samples)] <- paste0("SV", 1:svseq_1$n.sv)
ddsTxi_2 <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Status + Age + Ethnicity + SV1 + SV2)

# Cuantificar los genes previos al filtrado
length(counts(ddsTxi_2)) 

###################
# 2. Pre-filtrado #
###################

# Eliminar genes con bajo conteo del nuevo objeto (se recomienda especificar el tamaño del grupo más pequeño como el número mínimo de muestras).
smallestGroupSize <- 9
keep <- rowSums(counts(ddsTxi_2) >= 10) >= smallestGroupSize
ddsTxi_2 <- ddsTxi_2[keep,]

# Genes tras filtrado
length(counts(ddsTxi_2))

# Indicar grupo de referencia
ddsTxi_2$Status <- factor(ddsTxi_2$Status, levels = c("Control", "Asthma"))

##############################################################################################
# 3. Evaluación de la calidad de los datos mediante agrupamiento y visualización de muestras #
##############################################################################################

# Heatmap de las distancias entre muestras
# 1. Calcular la distancia entre muestras utilizando la transformación por varianza estabilizada (VST)
vsd <- vst(ddsTxi_2, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

# 2. Ajustar formato de los datos
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Status, vsd$Ethnicity, sep="-")
colnames(sampleDistMatrix) <- NULL

# 3.Generar una paleta de colores de tonos azules
colors <- colorRampPalette(c("#F0F8FF", "#0000FF", "#00008B", "#000080"))(255)

# 4.Generar heat map  
png("data/processed/Deseq2/women/sample_distance_women_heatmap.png", width = 800, height = 600)  
pheatmap(sampleDistMatrix, 
         clustering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)
dev.off()

# Gráfica de componentes principales de las muestras
# 1. Crear datos de pca 
pcaData <- plotPCA(vsd, intgroup=c("Status", "Ethnicity"), returnData=TRUE)

# 2. Calcular el porcentaje de varianza explicada por los dos primeros componentes principales
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 3. Generar gráfica de pca 
png("data/processed/Deseq2/women/pca_plot_women_per_sample.png", width = 800, height = 600)  
ggplot(pcaData, aes(PC1, PC2, color=Status, shape=Ethnicity)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()
dev.off()

# Emplear OUTRIDER (OUTlier in RNA-Seq fInDER)
#Permite detectar el número de genes con expresión aberrante por muestra. Si una muestra presenta un número muy alta debe valorarse su eliminación 

# 1. Extraer datos
countData<-counts(ddsTxi)
ods <- OutriderDataSet(countData=countData)

# Filtrar genes no expresados
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

# Correr pipeline de OUTRIDER entero
ods <- OUTRIDER(ods)

# Extraer resultados
res_aberrant <- results(ods)

# Obtener el número de genes aberrantes por muestra
aberrant_counts <- aberrant(ods, by="sample")

# Convertir a data frame para ggplot
df_aberrant <- data.frame(
  Sample = names(aberrant_counts),
  Aberrations = aberrant_counts
)

# Ordenar por número de aberraciones
df_aberrant <- df_aberrant[order(df_aberrant$Aberrations, decreasing = TRUE), ]

# Graficar el número de genes aberrantes por muestra
png("data/processed/Deseq2/women/aberrant_genes_women_per_sample.png", width = 800, height = 600)
ggplot(df_aberrant, aes(x = reorder(Sample, -Aberrations), y = Aberrations)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Número de genes aberrantes por muestra",
    x = "Muestras",
    y = "Número de genes aberrantes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#######################################
# 4. Análsis de expresión diferencial #
#######################################
# Realizar análisis con Deseq2
dds_2 <- DESeq(ddsTxi_2) # Por defecto se utiliza la prueba de Wald; también puede usarse la prueba de razón de verosimilitud (LRT) si se especifica. LRT se emplea para evaluar el efecto global de una covariable.

# Graficar la dispersión de los datos
png("data/processed/Deseq2/women/postfiltering_plotDispEsts_women_rnaseq_with_IMC.png", width = 800, height = 600)
plotDispEsts(dds_2)
dev.off()

# Ajustar resultados con ponderación de hipótesis independientes usando IHW
res <- results(dds_2, contrast=c("Status", "Asthma", "Control"), alpha=0.05, filterFun=ihw) 

# Anotar el nombre de los genes (symbol y entrez.id)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# Para ver los resultados usar:
# summary(res)

# Filtrar genes con p-valores ajustados significativos (padj < 0.05) y diferencialmente expresados (log2FoldChange > 0.05 o log2FoldChange < -0.05)
significativos <- res[(res$padj < 0.05 & (res$log2FoldChange > 0.05 | res$log2FoldChange < -0.05)), ]

# Guardar resultados significativos
write.table(significativos, "data/processed/Deseq2/women/DE_genes_0.05_women.txt", quote=F, row.names=T)

####################################
# 5. Exploración de los resultados #
####################################
# Para extraer variables usadas e información del test usar:
# mcols(res)$description

### OPCIONAL: Control de calidad de los datos de conteo a nivel de muestra si se reportan muchos valores atípicos (por ejemplo, cientos o miles) en el resumen (res), se podría considerar realizar una exploración adicional para ver si una muestra o algunas pocas muestras deben ser eliminadas debido a su baja calidad. REALIZADO POR COMPROBACIÓN EXTRA NO POR EXCESO DE VALORES ATÍPICOS.

# 1. Visualizar las distancias de Cook para la detección de valores atípicos en las muestras
#Verificar si las distancias de Cook de una muestra son consistentemente más altas que las de las demás.

par(mar=c(8,5,2,2))
png("data/processed/Deseq2/women/cooks_distance_sample_outlier_detection.png", width = 800, height = 600)
boxplot(log10(assays(dds_2)[["cooks"]]), range=0, las=2)
dev.off()

# 5.1 MA-plot
##############
# Extraer los valores desde el objeto de resultados de DESeq2
ma_data <- as.data.frame(res)

# Añadir una columna para resaltar diferencialmente expresados
ma_data$over_infra <- ma_data$log2FoldChange > 0.5 | ma_data$log2FoldChange < -0.5

# Graficar el MA plot 
png("data/processed/Deseq2/women/plotMA_rnaseq_women.png", width = 800, height = 600)
ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange, color = over_infra)) +
  geom_point(alpha = 0.5) +  # Puntos semi-transparentes
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Línea en log2FC = 0
  scale_color_manual(values = c("black", "blue")) +  # No significativo: negro, significativo: azul
  labs(
    title = "MA Plot (Manual - DESeq2)",
    x = "Log10 Base Mean",
    y = "Log2 Fold Change"
  ) +
  theme_minimal()
dev.off()

# 5.2 Volcano plot
####################
# Basado en el pipeline detallado en: https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

labels=res$symbol
png("data/processed/Deseq2/women/enhanced_volcano_plot_rnaseq_women.png", width = 800, height = 600)
EnhancedVolcano(res,
    lab = labels,
    x = 'log2FoldChange',
    y = 'pvalue',
title = 'Women volcano plot: Control vs Asthma',
pCutoff = 10e-5,
    FCcutoff = 0.5)
dev.off()

# 5.3 Análisis de enriquecimiento KEGG
#######################################
# Crear una lista con los nombres de los genes en Código Entrez.
genes <- significativos_01$entrez

# Crear una lista para los valores de expresión (FOLD CHANGE).
fold_change <- significativos_01$log2FoldChange

#Asignar los nombres de genes para cada resultado de expresión
names(fold_change)<-genes

#Generar un enriquecimiento con datos a partir de información del KEGG
KEGG_genes <- enrichKEGG(gene= genes, organism= "hsa", pvalueCutoff = 0.05)

# Gráfico de puntos
png("data/processed/Deseq2/women/KEGG_analysis_women_0.05.png", width = 800, height = 600)
dotplot(KEGG_genes, showCategory=20)
dev.off()

# Crear Emmaplot de enriquecimiento KEGG
KEGG_genes_sim <- pairwise_termsim(KEGG_genes)
png("data/processed/Deseq2/women/emapplot_routes_enrichment_analysis_women_0.05.png", width = 800, height = 600)
emapplot(KEGG_genes_sim, showCategory = 20)
dev.off()


# 5.4 Análisis de ontología génica (GO)
########################################
# Extraer datos 
significativos$ENSEMBL<-rownames(significativos)
go_data<-data.frame(significativos)
universe <- go_data$ENSEMBLE
sigGenes <- go_data %>% 
  filter(!is.na(ENSEMBL)) 
sigGenes<-sigGenes$entrez

# Análisis de enriquecimiento de procesos biológicos
enrich_go <- enrichGO(
  gene= sigGenes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  universe = universe,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable=TRUE
)

# Crear gráfica con los resultados significativos
png("data/processed/Deseq2/women/GO_analysis_women_0.05.png", width = 800, height = 600)
dotplot(enrich_go, showCategory = 20)
dev.off()


# Crear Emmaplot de enriquecimiento GO
enrich_go_sim <- pairwise_termsim(enrich_go)
png("data/processed/Deseq2/women/GO_emmaplot_analysis_women_0.05.png", width = 800, height = 600)
emapplot(enrich_go_sim, showCategory = 20)
dev.off()
```
