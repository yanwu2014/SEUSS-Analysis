library(methods)
library(Seurat)
library(perturbLM)

options <- commandArgs(trailingOnly = T)

## Command line parameters (can hardcode in script if desired)
matrix.file <- options[[1]] ## Input counts matrix (either in tab separated format or 10X genomics sparse format)
genotypes.file <- options[[2]] ## Genotype dictionary as a ".csv" file
output.file <- options[[3]] ## Output RData file
regress.batch <- as.logical(options[[4]]) ## Whether to run batch regression

## Other parameters
batch.exclude <- NULL ## Exclude a specific batch from analysis

min.cells.frac <- 0.01 ## Fraction of cells expressing gene in order to include
min.genes.exp <- 300 ## Minimum number of expressed genes for a cell
trim <- 0.005 ## Trim fraction

regress.model <- "negbinom" ## Model for regressing out unwanted effects
min.expr <- 0.1 ## Minimum average expression for variable genes
min.var.z <- 0.5 ## Minimum Z-score for variable genes
scale.var <- T ## Only scale variable genes

pcs.use <- 30 ## Number of PCs to use for clustering
cluster.resolution <- 1.8 ## Cluster resolution
min.conn <- 2e-4 ## Minimum SNN connectivity between clusters to calculate accuracy
acc.cutoff <- 0.9 ## Minimum classification accuracy to keep clusters separate
ntc <- "mCherry" ## Name of negative control genotype

# Load the dataset, filter and trim counts matrix
counts <- ReadData(matrix.file)
counts <- FilterData(counts, min.cells.frac = min.cells.frac, trim = trim, min.genes = min.genes.exp, min.expr = 0)

# Check for batch effects
if (regress.batch) {
  batch <- sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]])
  batch <- batch[!batch %in% batch.exclude]

  counts <- counts[,names(batch)]
  batch <- factor(batch)

  nUMI <- colSums(counts)
  pd <- data.frame(nUMI, batch)
} else {
  nUMI <- colSums(counts)
  pd <- data.frame(nUMI)
}
scale.factor <- median(nUMI)
print(paste("Scale factor:", scale.factor))

# Mitochondrial fraction
mito.genes <- grep("^MT-", rownames(counts), value = T)
percent.mt <- colSums(counts[mito.genes, ])/colSums(counts)
pd$percent.mt <- percent.mt
counts <- counts[!rownames(counts) %in% mito.genes,]

# Load genotypes
genotypes.list <- ReadGenotypes(genotypes.file)
genotypes.list <- lapply(genotypes.list, function(cells) cells[cells %in% colnames(counts)])
genotypes.list <- genotypes.list[sapply(genotypes.list, length) > 0]

guides.df <- as.matrix(DesignMatrixGenotypes(genotypes.list, max.genotypes = 3, min.cells = 25))
guides.df <- PadDesignMatrix(guides.df, colnames(counts))
guides.df <- data.frame(apply(guides.df, 2, function(x) factor(x)))
pd <- cbind(pd, guides.df)

# Run clustering pipeline via Seurat
se.obj <- CreateSeuratObject(raw.data = counts, project = "10x", normalization.method = "LogNormalize",
                             scale.factor = scale.factor, meta.data = pd)
print(dim(se.obj@data))

# Find overdispersed genes for PCA
se.obj <- FindVariableGenes(se.obj, x.low.cutoff = min.expr, x.high.cutoff = 8,
                            y.cutoff = min.var.z, y.high.cutoff = Inf, do.plot = T)

if (scale.var) {
  scale.genes <- se.obj@var.genes
} else {
  scale.genes <- rownames(se.obj@data)
}

# Regress out confounders
print(paste("Regress out batch effects:", regress.batch))
if (regress.batch) {
  se.obj <- ScaleData(se.obj, vars.to.regress = c("nUMI", "percent.mt", "batch"), model.use = regress.model,
                      scale.max = 10, genes.use = scale.genes)
} else {
  se.obj <- ScaleData(se.obj, vars.to.regress = c("nUMI", "percent.mt"), model.use = regress.model,
                      scale.max = 10, genes.use = scale.genes)
}

print("Done scaling data")
save.image(output.file)

se.obj <- RunPCA(se.obj, pc.genes = se.obj@var.genes, pcs.compute = 40, do.print = F)
PCElbowPlot(se.obj, num.pc = 40)
se.obj <- FindClusters(se.obj, reduction.type = "pca", dims.use = 1:pcs.use, resolution = cluster.resolution,
                       k.param = 30, save.SNN = T, prune.SNN = 1/20)

rm(counts)
se.obj@raw.data <- NULL

conn <- Seurat:::CalcConnectivity(se.obj)
print(sum(conn > min.conn))

se.obj <- ValidateClusters(se.obj, pc.use = 1:pcs.use, top.genes = 15, min.connectivity = min.conn, acc.cutoff = acc.cutoff, verbose = T)
se.obj <- BuildClusterTree(se.obj, do.reorder = T, reorder.numeric = T)

print("Done with clustering")
save.image(output.file)

se.obj <- RunTSNE(se.obj, reduction.use = "pca", do.fast = T, dims.use = 1:pcs.use)
TSNEPlot(se.obj)
if (regress.batch) {
  TSNEPlot(se.obj, group.by = "batch")
}

cluster.markers <- FindAllMarkers(se.obj, thresh.use = 0.1, test.use = "bimod", only.pos = T, print.bar = T)
cluster.markers <- subset(cluster.markers, p_val < 1e-2)

print("Done finding markers")
save.image(output.file)

# Cluster enrichment by genotype
clusters <- sapply(as.character(se.obj@ident), function(i) paste("C", i , sep = ""))
names(clusters) <- se.obj@cell.names

genotypes.list <- lapply(genotypes.list, function(x) x[x %in% se.obj@cell.names])
clusters.list <- UnflattenCellGenotypes(clusters)

df <- GenotypeClusterCounts(genotypes.list, clusters.list)
df.pvals <- GenotypeClusterPvals(df)

print("Done with genotype enrichment")
save.image(output.file)

saveRDS(se.obj, file = gsub(".RData", ".Robj", output.file))
write.table(df, file = gsub("clustering.RData", "cluster_genotype_counts.tsv", output.file), sep = "\t")
write.table(df.pvals, file = gsub("clustering.RData", "cluster_enrichment_pvals.tsv", output.file), sep = "\t")
