if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Matriz de counts
counts <- read.csv("data/counts_matrix.csv", row.names = 1)

# Matriz de metadata
meta <- read.csv("data/sample_metadata.csv", row.names = 1)

# Comprobamos el orden 
colnames(counts)
rownames(meta)

all(colnames(counts)) == rownames(meta))

meta$condition <- factor(meta$condition, levels = c("control", "LPS"))
meta

# Crear el objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ condition
)

# Ejecutar el pipeline interno de DESeq2
dds <- DESeq(dds)

# Obtener resultados DE entre LPS y control
res <- results(dss)
head(res)

# Ordenar por p-valor
res_ordenado <- res[order(res$pvalue), ]

# Ver los genes mas significativos
head(res_ordenado)

# Salida CSV
res_df <- as.data.frame(res_ordered)

write.csv(
  res_df,
  file = "results/deseq2_results_microglia_simple.csv"
)




