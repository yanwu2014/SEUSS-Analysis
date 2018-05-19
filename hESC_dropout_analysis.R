library(matrixStats)
library(swne)
library(perturbLM)
library(ggplot2)

## Load Data
counts <- read.table("up-tf_dropout_counts.tsv", sep = "\t", header = T)
rownames(counts) <- counts[[1]]
counts <- as.matrix(counts[,3:ncol(counts)])

## Normalize counts
norm.counts <- counts/colSums(counts)

## Calculate log-fold change against control
log.fc.reads <- log2(apply(norm.counts, 2, function(x) x/norm.counts[,"ctrl"]))
log.fc.reads <- log.fc.reads[,!colnames(log.fc.reads) == "ctrl"]

## Dropout analysis with genome reads
heat.mat.reads <- as.matrix(data.frame(Stem = rowMeans(log.fc.reads[,c("high1", "high2", "low1a", "low1b")]), 
                                       Endo = log.fc.reads[,"egm"], Multi = log.fc.reads[,"dmem"]))
heat.mat.reads <- heat.mat.reads[order(rowMeans(heat.mat.reads), decreasing = F),]
colnames(heat.mat.reads) <- c("Pluripotent", "Endothelial", "Multilineage")

pdf("up-tf_dropout_reads_heatmap.pdf", width = 2, height = 8)
ggHeat(heat.mat.reads, clustering = "col", x.lab.size = 12, y.lab.size = 6,
       heatscale = c(low = 'purple', mid = 'white', high = 'green')) + scale_color_brewer(direction = -1)
dev.off()


## Dropout analysis with cell counts
genotype.files <- c("Individual_Batch_Analysis/up-tf2_pheno_dict.csv",
                    "Individual_Batch_Analysis/up-tf3_pheno_dict.csv",
                    "Individual_Batch_Analysis/up-tf8_pheno_dict.csv",
                    "Individual_Batch_Analysis/up-tf10_pheno_dict.csv",
                    "Individual_Batch_Analysis/up-tf11_pheno_dict.csv",
                    "Individual_Batch_Analysis/up-tf12_pheno_dict.csv")

cell.counts.list <- lapply(genotype.files, function(x) sapply(ReadGenotypes(x), length))
common.genotypes <- Reduce(intersect, lapply(cell.counts.list, names))

## Create cell counts matrix
cell.counts <- do.call("cbind", lapply(cell.counts.list, function(x) x[common.genotypes]))
colnames(cell.counts) <- c("stem1", "stem2", "egm1", "dmem1", "egm2", "dmem2")

## Normalize cell counts
norm.cell.counts <- cell.counts/colSums(cell.counts)

## Calculate logFC against plasmid control
ctrl.counts <- norm.counts[rownames(cell.counts), "ctrl"]

log.fc.cells <- log2(apply(norm.cell.counts, 2, function(x) x/ctrl.counts))
heat.mat.cells <- as.matrix(data.frame("Stem" = rowMeans(log.fc.cells[,c("stem1","stem2")]),
                                       "Endo" = rowMeans(log.fc.cells[,c("egm1","egm2")]),
                                       "Multi" = rowMeans(log.fc.cells[,c("dmem1","dmem2")])))
heat.mat.cells <- heat.mat.cells[order(rowMeans(heat.mat.cells), decreasing = F),]
colnames(heat.mat.cells) <- c("Pluripotent", "Endothelial", "Multilineage")
heat.mat.cells <- heat.mat.cells[rownames(heat.mat.cells) != "OTX2",]

pdf("up-tf_dropout_cells_heatmap.pdf", width = 2, height = 8)
ggHeat(heat.mat.cells, clustering = "col", x.lab.size = 12, y.lab.size = 6,
       heatscale = c(low = 'purple', mid = 'white', high = 'green')) + scale_color_brewer(direction = -1)
dev.off()


## Correlate fitness replicates from read counts
genes.lbl <- abs(log.fc.reads[,"high1"]) > 2.5 | abs(log.fc.reads[,"high2"]) > 2.5
pdf("up-tf-stem_fitness_reads_replicate_corr.pdf", width = 4.5, height = 4.5)
PlotCorrelation(log.fc.reads[,"high1"], log.fc.reads[,"high2"],
                "Replicate 1", "Replicate 2", "Pluripotent gDNA logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()


## Correlate fitness replicates from cell counts
log.fc.cells <- log.fc.cells[rownames(log.fc.cells) != "OTX2",]
genes.lbl <- abs(log.fc.cells[,"stem1"]) > 2.5 | abs(log.fc.cells[,"stem2"]) > 2.5
pdf("up-tf-stem_fitness_cells_replicate_corr.pdf", width = 4.5, height = 4.5)
PlotCorrelation(log.fc.cells[,"stem1"], log.fc.cells[,"stem2"],
                "Replicate 1", "Replicate 2", "Pluripotent scRNA logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()

genes.lbl <- abs(log.fc.cells[,"egm1"]) > 2.75 | abs(log.fc.cells[,"egm2"]) > 2.75
pdf("up-tf-endo_fitness_cells_replicate_corr.pdf", width = 4.5, height = 4.5)
PlotCorrelation(log.fc.cells[,"egm1"], log.fc.cells[,"egm2"],
                "Replicate 1", "Replicate 2", "Endothelial scRNA logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()

genes.lbl <- abs(log.fc.cells[,"dmem1"]) > 4 | abs(log.fc.cells[,"dmem2"]) > 4
pdf("up-tf-multi_fitness_cells_replicate_corr.pdf", width = 4.5, height = 4.5)
PlotCorrelation(log.fc.cells[,"dmem1"], log.fc.cells[,"dmem2"],
                "Replicate 1", "Replicate 2", "Multilineage scRNA logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()


## Correlate fitness from genome reads with fitness from cell counts
pdf("up-tf-stem_cells_vs_reads_dropout.pdf", height = 4, width = 4.5)
genes.lbl <- abs(log.fc.cells[,"Pluripotent"]) > 2
PlotCorrelation(log.fc.cells[,"Pluripotent"], log.fc.reads[rownames(log.fc.cells), "Pluripotent"], 
                "scRNA logFC", "gDNA logFC", "Fitness Correlation", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()

pdf("up-tf-endo_cells_vs_reads_dropout.pdf", height = 4, width = 4.5)
genes.lbl <- abs(log.fc.cells[,"Endothelial"]) > 1.95
PlotCorrelation(log.fc.cells[,"Endothelial"], log.fc.reads[rownames(log.fc.cells), "Endothelial"], 
                "scRNA logFC", "gDNA logFC", "Fitness Correlation", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()

pdf("up-tf-multi_cells_vs_reads_dropout.pdf", height = 4, width = 4.5)
genes.lbl <- abs(log.fc.cells[,"Multilineage"]) > 3.5
PlotCorrelation(log.fc.cells[,"Multilineage"], log.fc.reads[rownames(log.fc.cells), "Multilineage"], 
                "scRNA logFC", "gDNA logFC", "Fitness Correlation", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T)
dev.off()


## Correlate differential genes with fitness 
# summary.file <- "up-tf-stem_summary.tsv"
# summary.file <- "up-tf-endo_summary.tsv"
summary.file <- "up-tf-multi_summary.tsv"

summary.df <- read.table(summary.file, sep = "\t", header = T, row.names = 1)
n.genes <- summary.df$n_diff_genes
names(n.genes) <- rownames(summary.df)

# log.fc.use <- log.fc.reads[,"Pluripotent"]
# log.fc.use <- log.fc.reads[,"Endothelial"]
log.fc.use <- log.fc.reads[,"Multilineage"]

genotypes.use <- intersect(names(n.genes), names(log.fc.use))
n.genes <- n.genes[genotypes.use]; log.fc.use <- log.fc.use[genotypes.use];

pts.label <- n.genes > 60 | log.fc.use > 1 | log.fc.use < -4
print(sum(pts.label))

pdf("up-tf-multi_diff_genes_vs_fitness.pdf", width = 4, height = 4)
PlotCorrelation(n.genes, log.fc.use, NULL, NULL, show.corr = F, box = F, use.label = T, pts.label = pts.label,
                title = "Multilineage media", pt.size = 3.0)
dev.off()