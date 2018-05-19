library(Seurat)
library(swne)
library(perturbLM)
library(cellMapper)

## Set parameters
file.handle <- "up-tf-stem"
# file.handle <- "up-tf-endo"
# file.handle <- "up-tf-multi"

min.cells.frac <- 0.001
min.genes.exp <- 200
trim <- 0.001
pcs.use <- 25
min.genotype.cells <- 20

## Load counts
counts <- ReadData(file.handle)

## Load genotypes
genotypes.list <- ReadGenotypes(paste(file.handle, "pheno_dict.csv", sep = "_"))

## Load genotype summary
min.cl.pval <- 1e-12
min.diff.genes <- 50

genotype.summary <- read.table(paste(file.handle, "summary.tsv", sep = "_"), sep = "\t", header = T, row.names = 1)
sig.genotypes <- union(rownames(subset(genotype.summary, abs(min_cluster_pval) < min.cl.pval)),
                       rownames(subset(genotype.summary, n_diff_genes >= min.diff.genes)))
sig.genotypes <- sort(unique(c(sig.genotypes, "mCherry-int-ctrl"))); print(sig.genotypes);

## Filter counts
counts <- counts[,colnames(counts) %in% unlist(genotypes.list, F, F)]

## Filter genotypes
genotypes.list <- lapply(genotypes.list, function(cells) cells[cells %in% colnames(counts)])
genotypes.list <- genotypes.list[sapply(genotypes.list, length) > 0]

guides.matrix <- DesignMatrixGenotypes(genotypes.list, max.genotypes = 2, min.cells = 5, drop.cells = T)
guides.matrix <- CleanDesignCtrl(guides.matrix, ctrl = "mCherry")

n_cells <- Matrix::colSums(guides.matrix); n_cells <- n_cells[n_cells >= min.genotype.cells];
guides.matrix <- guides.matrix[,names(n_cells)]
guides.matrix <- guides.matrix[Matrix::rowSums(guides.matrix) > 0,]

## Filter and scale counts
counts <- counts[,rownames(guides.matrix)]
counts <- FilterData(counts, min.cells.frac, trim, min.genes.exp)
dim(counts)

batch <- factor(sapply(colnames(counts), ExtractField, field = 2, delim = "\\."))
norm.counts <- ScaleCounts(counts, batch = batch, method = "log", adj.var = T)

## Get TF marker genes
max.fdr <- 0.05
min.abs.cf <- 0.025

coefs.df <- read.table(paste(file.handle, "regression.pvals.tsv", sep = "_"), sep = "\t", header = T)
coefs.df <- subset(coefs.df, abs(cf) > min.abs.cf & FDR < max.fdr)

## Load reference data
ref.obj <- readRDS("Reference_Data/mouse-cell-atlas_seurat.Robj")

## Select genes to use for classification
ref.markers.df <- read.table("Reference_Data/mouse-cell-atlas_markers.txt", sep = "\t", header = T)
mouse2human <- MouseHumanMapping(as.character(ref.markers.df$gene))
ref.markers.df <- subset(ref.markers.df, gene %in% names(mouse2human))
ref.markers.df$gene <- mouse2human[ref.markers.df$gene]

genes.use <- union(as.character(ref.markers.df$gene), unique(coefs.df$Gene))
genes.use <- Reduce(intersect, list(genes.use, rownames(ref.obj@data), rownames(norm.counts)))
length(genes.use)

## Scale reference and teratoma data
ref.norm.counts <- ScaleCounts(ref.obj@raw.data[,ref.obj@cell.names], method = "log")

## Run PCA on reference data
ref.pca <- FastPCA(ref.norm.counts, genes.use, pcs.use)
ref.pc.emb <- ref.pca$pc.emb

## Project teratoma onto reference PCs
pc.emb <- ProjectPCA(norm.counts, ref.pca$pc.load, genes.use)

## Classify cells with kNN classifier
ref.clusters <- ref.obj@meta.data$tissue; names(ref.clusters) <- ref.obj@cell.names; levels(ref.clusters)
mapped.clusters <- MapCellsKNN(pc.emb, ref.pc.emb, ref.ident = ref.clusters, k = 40)

## Summarize cell classifications by genotype
genotypes.list <- lapply(sig.genotypes, function(g) rownames(guides.matrix[guides.matrix[,g] == 1,]))
names(genotypes.list) <- sig.genotypes
genotype.frac <- SummarizeCellMapping(mapped.clusters, genotypes.list, ident.cell.thresh = 0.5,
                                      hard.cell.thresh = T, min.ref.frac = 0.1)

pdf(paste(file.handle, "mca_mapping.pdf", sep = "_"), width = 5.5, height = 3.5)
ggHeat(genotype.frac, clustering = "both")
dev.off()
