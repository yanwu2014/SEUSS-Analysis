library(methods)
library(matrixStats)
library(perturbLM)

## Load Data
raw.counts <- read.table("up-tf_dropout_counts.tsv", sep = "\t", header = T)
rownames(raw.counts) <- raw.counts[[1]]; raw.counts <- as.matrix(raw.counts[,3:ncol(raw.counts)]);
counts <- raw.counts/colSums(raw.counts)

log.fc <- log2(apply(counts, 2, function(x) x/counts[,"ctrl"]))
log.fc <- log.fc[,!colnames(log.fc) == "ctrl"]

## Summary figure with genome reads
log.fc.mtesr <- rowMeans(log.fc[,c("high1", "high2", "low1a", "low1b")])
n_reads_logfc <- as.matrix(data.frame(Stem = log.fc.mtesr, Endo = log.fc[,"egm"], Multi = log.fc[,"dmem"]))
n_reads_logfc <- n_reads_logfc[order(rowMeans(n_reads_logfc), decreasing = F),]
colnames(n_reads_logfc) <- c("Stem cell", "Endothelial", "Multilineage")

# pdf("dropout_heatmap.pdf", width = 2, height = 8)
ggheat(n_reads_logfc, clustering = "col", x.lab.size = 12, y.lab.size = 6,
       heatscale = c(low = 'purple', mid = 'white', high = 'green')) + scale_color_brewer(direction = -1)
# dev.off()


## Summary figure with cell counts
genotype_files <- c("hESC_Genotype_Dictionaries/up-tf2_pheno_dict.csv",
                    "hESC_Genotype_Dictionaries/up-tf3_pheno_dict.csv",
                    "hESC_Genotype_Dictionaries/up-tf8_pheno_dict.csv",
                    "hESC_Genotype_Dictionaries/up-tf10_pheno_dict.csv",
                    "hESC_Genotype_Dictionaries/up-tf11_pheno_dict.csv",
                    "hESC_Genotype_Dictionaries/up-tf12_pheno_dict.csv")

n_cells_counts_list <- lapply(genotype_files, function(x) sapply(ReadGenotypes(x), length))
common_tfs <- Reduce(intersect, lapply(n_cells_counts_list, function(x) names(x)))
n_cells_counts <- do.call("cbind", lapply(n_cells_counts_list, function(x) x[common_tfs]))
colnames(n_cells_counts) <- c("stem1", "stem2", "egm1", "dmem1", "egm2", "dmem2")
n_cells_counts <- n_cells_counts/colSums(n_cells_counts)

ctrl_counts <- counts[rownames(n_cells_counts), "ctrl"]
n_cells_logfc <- log2(apply(n_cells_counts, 2, function(x) x/ctrl_counts))
n_cells_logfc <- as.matrix(data.frame("Stem" = rowMeans(n_cells_logfc[,c("stem1","stem2")]),
                                      "Endo" = rowMeans(n_cells_logfc[,c("egm1","egm2")]),
                                      "Multi" = rowMeans(n_cells_logfc[,c("dmem1","dmem2")])))
n_cells_logfc <- n_cells_logfc[order(rowMeans(n_cells_logfc), decreasing = F),]
colnames(n_cells_logfc) <- c("Stem cell", "Endothelial", "Multilineage")
n_cells_logfc <- n_cells_logfc[rownames(n_cells_logfc) != "OTX2",]

# pdf("cell_counts_dropout_heatmap.pdf", width = 2, height = 8)
ggheat(n_cells_logfc, clustering = "col", x.lab.size = 12, y.lab.size = 6,
       heatscale = c(low = 'purple', mid = 'white', high = 'green')) + scale_color_brewer(direction = -1)
# dev.off()


## Correlate fitness from genome reads with fitness from cell counts
# pdf("hPSC_media_cell_counts_vs_dropout.pdf", height = 4, width = 4.5)
genes.lbl <- abs(n_cells_logfc[,"Stem cell"]) > 2
ggcorrelation(n_cells_logfc[,"Stem cell"], n_reads_logfc[rownames(n_cells_logfc), "Stem cell"], 
              "scRNA fitness logFC", "gDNA fitness logFC", "scRNA vs gDNA", pt.size = 1.5,
              pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
# dev.off()

# pdf("endothelial_media_cell_counts_vs_dropout.pdf", height = 4, width = 4.5)
genes.lbl <- abs(n_cells_logfc[,"Endothelial"]) > 1.95
ggcorrelation(n_cells_logfc[,"Endothelial"], n_reads_logfc[rownames(n_cells_logfc), "Endothelial"], 
              "scRNA fitness logFC", "gDNA fitness logFC", "scRNA vs gDNA", pt.size = 1.5,
              pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
# dev.off()

# pdf("multilineage_media_cell_counts_vs_dropout.pdf", height = 4, width = 4.5)
genes.lbl <- abs(n_cells_logfc[,"Multilineage"]) > 3.5
ggcorrelation(n_cells_logfc[,"Multilineage"], n_reads_logfc[rownames(n_cells_logfc), "Multilineage"], 
              "scRNA fitness logFC", "gDNA fitness logFC", "scRNA vs gDNA", pt.size = 1.5,
              pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
# dev.off()


## Correlate differential genes with fitness 
# summary.file <- "up-tf-all_genotype_summary.tsv"
# summary.file <- "up-tf-endo_genotype_summary.tsv"
summary.file <- "up-tf-multi_genotype_summary.tsv"
summary.df <- read.table(summary.file, sep = "\t")

n_diff_genes <- summary.df$n_diff_genes; names(n_diff_genes) <- rownames(summary.df);
# log.fc.use <- n_reads_logfc[,"Stem cell"]
# log.fc.use <- n_reads_logfc[,"Endothelial"]
log.fc.use <- n_reads_logfc[,"Multilineage"]

genotypes.use <- intersect(names(n_diff_genes), names(log.fc.use))
n_diff_genes <- n_diff_genes[genotypes.use]; log.fc.use <- log.fc.use[genotypes.use];

top.pts <- c(names(n_diff_genes[n_diff_genes > 500]), names(log.fc.use[abs(log.fc.use) > 4]))
pts.label <- names(n_diff_genes) %in% top.pts
print(sum(pts.label))

# pdf("multi_diff_genes_vs_fitness.pdf", width = 3.5, height = 3.5)
ggcorrelation(n_diff_genes, log.fc.use, NULL, NULL, show.corr = F, box = F, use.label = T, pts.label = pts.label,
              title = "Multilineage media")
# dev.off()