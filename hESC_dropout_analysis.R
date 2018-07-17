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
       heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato')) + scale_color_brewer(direction = -1)
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
cell.counts <- cell.counts[rownames(cell.counts) != "OTX2",]

heat.mat.cells.counts <- as.matrix(data.frame("Stem" = rowSums(cell.counts[,c("stem1","stem2")]),
                                              "Endo" = rowSums(cell.counts[,c("egm1","egm2")]),
                                              "Multi" = rowSums(cell.counts[,c("dmem1","dmem2")])))

## Normalize cell counts
norm.cell.counts <- cell.counts/colSums(cell.counts)
heat.mat.cells.counts <- heat.mat.cells.counts/colSums(heat.mat.cells.counts)

## Calculate logFC against plasmid control
ctrl.counts <- norm.counts[rownames(cell.counts), "ctrl"]

log.fc.cells <- log2(apply(norm.cell.counts, 2, function(x) x/ctrl.counts))
heat.mat.cells <- log2(apply(heat.mat.cells.counts, 2, function(x) x/ctrl.counts))
heat.mat.cells <- heat.mat.cells[order(rowMeans(heat.mat.cells), decreasing = F),]
colnames(heat.mat.cells) <- c("Pluripotent", "Endothelial", "Multilineage")

pdf("up-tf_dropout_cells_heatmap.pdf", width = 2, height = 8)
ggHeat(heat.mat.cells, clustering = "col", x.lab.size = 12, y.lab.size = 6,
       heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato')) + scale_color_brewer(direction = -1)
dev.off()


## Correlate fitness replicates from read counts
genes.lbl <- abs(log.fc.reads[,"high1"]) > 2.5 | abs(log.fc.reads[,"high2"]) > 2.5
pdf("up-tf-stem_fitness_reads_replicate_corr.pdf", width = 3.25, height = 3.25)
PlotCorrelation(log.fc.reads[,"high1"], log.fc.reads[,"high2"],
                "Replicate 1", "Replicate 2", "Pluripotent gDNA logFC", pt.size = 2.0, 
                pts.label = genes.lbl, use.label = T, 
                box = T, font.size = 12, show.corr = T)
dev.off()


## Correlate fitness replicates from cell counts
genes.lbl <- abs(log.fc.cells[,"stem1"]) > 3 | abs(log.fc.cells[,"stem2"]) > 3
pdf("up-tf-stem_fitness_cells_replicate_corr.pdf", width = 3.25, height = 3.25)
PlotCorrelation(log.fc.cells[,"stem1"], log.fc.cells[,"stem2"],
                "Replicate 1", "Replicate 2", "logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T,
                font.size = 12)
dev.off()

genes.lbl <- abs(log.fc.cells[,"egm1"]) > 3 | abs(log.fc.cells[,"egm2"]) > 3
pdf("up-tf-endo_fitness_cells_replicate_corr.pdf", width = 3.25, height = 3.25)
PlotCorrelation(log.fc.cells[,"egm1"], log.fc.cells[,"egm2"],
                "Replicate 1", "Replicate 2", "logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T,
                font.size = 12)
dev.off()

genes.lbl <- abs(log.fc.cells[,"dmem1"]) > 4.5 | abs(log.fc.cells[,"dmem2"]) > 4.5
pdf("up-tf-multi_fitness_cells_replicate_corr.pdf", width = 3.25, height = 3.25)
PlotCorrelation(log.fc.cells[,"dmem1"], log.fc.cells[,"dmem2"],
                "Replicate 1", "Replicate 2", "logFC", pt.size = 3.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T,
                font.size = 12)
dev.off()


## Correlate fitness from genome reads with fitness from cell counts
pdf("up-tf-stem_cells_vs_reads_dropout.pdf", height = 3.25, width = 3.25)
genes.lbl <- abs(heat.mat.cells[,"Pluripotent"]) > 2.5
PlotCorrelation(heat.mat.cells[,"Pluripotent"], heat.mat.reads[rownames(heat.mat.cells), "Pluripotent"], 
                "scRNA logFC", "gDNA logFC", "Pluripotent fitness", pt.size = 2.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T, font.size = 12)
dev.off()

pdf("up-tf-endo_cells_vs_reads_dropout.pdf", height = 3.25, width = 3.25)
genes.lbl <- abs(heat.mat.cells[,"Endothelial"]) > 1.95
PlotCorrelation(heat.mat.cells[,"Endothelial"], heat.mat.reads[rownames(heat.mat.cells), "Endothelial"], 
                "scRNA logFC", "gDNA logFC", "Endothelial Fitness", pt.size = 2.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T, font.size = 12)
dev.off()

pdf("up-tf-multi_cells_vs_reads_dropout.pdf", height = 3.25, width = 3.25)
genes.lbl <- abs(heat.mat.cells[,"Multilineage"]) > 3
PlotCorrelation(heat.mat.cells[,"Multilineage"], heat.mat.reads[rownames(heat.mat.cells), "Multilineage"], 
                "scRNA logFC", "gDNA logFC", "Multilineage Fitness", pt.size = 2.0,
                pts.label = genes.lbl, use.label = T, box = T, show.corr = T, font.size = 12)
dev.off()


## Correlate differential genes with fitness 
summary.file <- "up-tf-stem_summary.tsv"
# summary.file <- "up-tf-endo_summary.tsv"
# summary.file <- "up-tf-multi_summary.tsv"

summary.df <- read.table(summary.file, sep = "\t", header = T, row.names = 1)
n.genes <- summary.df$n_diff_genes
names(n.genes) <- rownames(summary.df)

log.fc.use <- heat.mat.reads[,"Pluripotent"]
# log.fc.use <- heat.mat.reads[,"Endothelial"]
# log.fc.use <- heat.mat.reads[,"Multilineage"]

genotypes.use <- intersect(names(n.genes), names(log.fc.use))
n.genes <- n.genes[genotypes.use]; log.fc.use <- log.fc.use[genotypes.use];

pts.label <- n.genes > 110 | log.fc.use > 2 | log.fc.use < -1.75
print(sum(pts.label))

pdf("up-tf-stem_diff_genes_vs_fitness.pdf", width = 4, height = 4)
PlotCorrelation(n.genes, log.fc.use, NULL, NULL, show.corr = F, box = F, use.label = T, pts.label = pts.label,
                title = "Pluripotent media", pt.size = 3.0)
dev.off()