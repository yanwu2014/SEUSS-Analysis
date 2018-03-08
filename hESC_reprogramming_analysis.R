library(swne)
library(Seurat)

counts <- ReadData("up-tf-all/")

genotypes.list <- ReadGroups("up-tf-all_pheno_dict.csv", sep.char = ",")
genotypes.list <- lapply(genotypes.list, function(x) x[x %in% colnames(counts)])
counts <- counts[,colnames(counts) %in% unlist(genotypes.list, F, F)]

tf.counts <- ReadData("up-tf-all_tf_barcodes.tsv", make.sparse = F)
tf.counts <- tf.counts[,colnames(counts)]

tfs.use <- c("KLF4", "CDX2", "SNAI2", "mCherry")

batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
head(batch)

genotypes.list <- genotypes.list[tfs.use]
genotypes <- FlattenGroups(genotypes.list)

counts.use <- FilterData(counts[,names(genotypes)], min.samples.frac = 0.025, trim = 0.001, 
                         min.nonzero.features = 200)
genotypes <- genotypes[colnames(counts.use)]
cells.use <- names(genotypes); batch <- batch[cells.use];

tf.counts.use <- tf.counts[,colnames(counts.use)]
for (tf in rownames(tf.counts.use)) {
  if (tf %in% rownames(counts.use)) {
    tf.counts.use[tf,] <- tf.counts.use[tf,] + counts.use[tf,]
  }
}
tf.counts.use <- tf.counts.use/colSums(tf.counts.use) * median(colSums(tf.counts.use))
tf.counts.use <- log(tf.counts.use + 1)
rownames(tf.counts.use) <- paste(rownames(tf.counts.use), "TF", sep = "_")

## Create Seurat object
meta.data <- data.frame(t(tf.counts.use))
meta.data$Batch <- batch
meta.data$Genotype <- genotypes

se.obj <- CreateSeuratObject(counts.use, meta.data = meta.data)
se.obj <- NormalizeData(se.obj, scale.factor = median(Matrix::colSums(counts.use)))
se.obj <- SetAllIdent(se.obj, "Genotype")

## Figure out which genes to use in the analysis
gene.module.clusters <- ReadGenesets("up-tf_functional_gene_modules.gmt")
gene.module.mapping <- read.table("up-tf_gene_module_mapping.txt", sep = "\t", header = T)

pluripotency.genes <- gene.module.clusters[["GM11"]]
se.obj <- FindVariableGenes(se.obj, x.low.cutoff = 0.05, y.cutoff = 1.5)
var.genes <- se.obj@var.genes
genes.use <- union(pluripotency.genes, var.genes)
# genes.use <- pluripotency.genes
genes.use <- genes.use[genes.use %in% rownames(counts.use)]
length(genes.use)

se.obj <- ScaleData(se.obj, genes.use = NULL, model.use = "negbinom", vars.to.regress = c("nUMI", "Batch"))
se.obj <- RunPCA(se.obj, genes.use = genes.use, pcs.compute = 30, do.print = F)
PCElbowPlot(se.obj, num.pc = 30)

## Monocle2 analysis
library(monocle)

order.genes <- genes.use

pd <- new("AnnotatedDataFrame", data = se.obj@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = order.genes))
rownames(fd) <- order.genes

## create a CDS with data.info.genes 
cds <- newCellDataSet(as.matrix(counts.use[order.genes,]), 
                      phenoData = pd, featureData = fd, lowerDetectionLimit = 0,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
plot_pc_variance_explained(cds)

## DDRTree reduction and cell ordering
cds@phenoData@data$batch <- batch[rownames(cds@phenoData@data)]
cds <- reduceDimension(cds, verbose = F, residualModelFormulaStr = ~batch, max_components = 10)
cds <- orderCells(cds, root_state = NULL, reverse = T)
save.image("up-tf-all_regrogramming.RData")

plot_cell_trajectory(cds, color_by = "State")
cds <- orderCells(cds, root_state = 15, reverse = T)

pdf("reprogramming_cell_trajectory.pdf", width = 6, height = 6)
plot_cell_trajectory(cds, color_by = "Genotype", show_branch_points = F)
dev.off()

pdf("reprogramming_pseudotime_overlay.pdf", width = 6, height = 6)
plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = F) + 
  scale_color_distiller(palette = "YlOrRd", direction = 1, 
                        guide = guide_colorbar(title = "Pseudotime", ticks = T, label = T))
dev.off()

## Pull out pseudotime
ps.time <- cds@phenoData@data$Pseudotime; names(ps.time) <- rownames(cds@phenoData);
se.obj <- AddMetaData(se.obj, ps.time, col.name = "Pseudotime")

pdf("reprogramming_pseudotime_vlnplot.pdf", width = 4, height = 3)
VlnPlot(se.obj, "Pseudotime", point.size.use = -1)
dev.off()

## Correlate genes with pseudotime
gene.ps.cors <- apply(se.obj@scale.data, 1, function(x) cor(ps.time, x))
tf.ps.cors <- apply(tf.counts.use, 1, function(x) cor(ps.time, x))
gene.ps.cors <- gene.ps.cors[!grepl("RPS|RPL|RP11|MT-|LINC", names(gene.ps.cors))]
gene.ps.cors <- c(tf.ps.cors, gene.ps.cors); gene.ps.cors <- sort(gene.ps.cors, decreasing = T)

head(gene.ps.cors, n = 50)
tail(gene.ps.cors, n = 50)

genes.plot <- c("DNMT3B", "POU5F1", "SOX2", "SALL2", "THY1", "TERF1",
                "LEFTY2", "LEFTY1", "EPCAM", "KLF4_TF", "CDX2_TF")
genes.cors <- sort(gene.ps.cors[genes.plot], decreasing = F)
barplot.df <- data.frame(R = genes.cors, Gene = factor(names(genes.cors), levels = names(genes.cors)))

pdf("reprogramming_pseudotime_barplot.pdf", width = 5.5, height = 4.0)
ggplot(data = barplot.df, aes(x = Gene, y = R)) + 
  geom_bar(fill = "skyblue4", position = "dodge", stat = "identity", width = 0.6) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) + ylab("Pearson")
dev.off()

save.image("up-tf-all_regrogramming.RData")

gene <- "CDX2_TF"
gene.expr <- tf.counts.use[gene,]

pdf(paste("reprogramming", gene, "vlnplot.pdf", sep = "_"), width = 4, height = 3)
VlnPlot(se.obj, gene, point.size.use = -1) + ggtitle("CDX2 Expr")
dev.off()


