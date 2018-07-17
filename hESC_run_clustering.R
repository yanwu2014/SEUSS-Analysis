library(methods)
library(Seurat)
library(swne)

options <- commandArgs(trailingOnly = T)

## Command line parameters (can hardcode in script if desired)
matrix.file <- options[[1]] ## Input counts matrix (either in tab separated format or 10X genomics sparse format)
output.file <- options[[2]] ## Output RData file
regress.batch <- as.logical(options[[3]]) ## Whether to run batch regression

## Other parameters
min.cells.frac <- 0.01 ## Fraction of cells expressing gene in order to include
min.genes.exp <- 300 ## Minimum number of expressed genes for a cell
trim <- 0.005 ## Trim fraction

regress.model <- "negbinom" ## Model for regressing out unwanted effects
min.expr <- 0.1 ## Minimum average expression for variable genes
min.var.z <- 0.5 ## Minimum Z-score for variable genes

pcs.use <- 30 ## Number of PCs to use for clustering
cluster.res <- 1.8 ## Cluster resolution
min.conn <- 1e-4 ## Minimum SNN connectivity between clusters to calculate accuracy
acc.cutoff <- 0.95 ## Minimum classification accuracy to keep clusters separate

## Load the dataset, filter and trim counts matrix
counts <- ReadData(matrix.file)
counts <- FilterData(counts, min.cells.frac,trim, min.genes.exp)
nUMI <- Matrix::colSums(counts)

## Check for batch effects
if (regress.batch) {
  batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
  pd <- data.frame(nUMI, batch)
} else {
  pd <- data.frame(nUMI)
}

## Mitochondrial fraction
mito.genes <- grep("^MT-", rownames(counts), value = T)
pd$percent.mt <- Matrix::colSums(counts[mito.genes, ])/Matrix::colSums(counts)
counts <- counts[!rownames(counts) %in% mito.genes,]

## Run clustering pipeline via Seurat
se.obj <- CreateSeuratObject(raw.data = counts, project = "10x", normalization.method = "LogNormalize",
                             scale.factor = median(nUMI), meta.data = pd)
dim(se.obj@data)

## Find overdispersed genes for PCA
se.obj <- FindVariableGenes(se.obj, x.low.cutoff = min.expr, x.high.cutoff = 8,
                            y.cutoff = min.var.z, y.high.cutoff = Inf, do.plot = T)
length(se.obj@var.genes)

## Regress out confounders
print(paste("Regress out batch effects:", regress.batch))
if (regress.batch) {
  se.obj <- ScaleData(se.obj, vars.to.regress = c("nUMI", "percent.mt", "batch"), model.use = regress.model,
                      scale.max = 10, genes.use = se.obj@var.genes)
} else {
  se.obj <- ScaleData(se.obj, vars.to.regress = c("nUMI", "percent.mt"), model.use = regress.model,
                      scale.max = 10, genes.use = se.obj@var.genes)
}

print("Done scaling data")
save.image(output.file)

## Run PCA and clustering
se.obj <- RunPCA(se.obj, pc.genes = se.obj@var.genes, pcs.compute = 40, do.print = F)
se.obj <- FindClusters(se.obj, dims.use = 1:pcs.use, resolution = cluster.res,
                       k.param = 30, save.SNN = T, prune.SNN = 1/20, print.output = F)

## Validate clustering
conn <- Seurat:::CalcConnectivity(se.obj)
print(sum(conn > min.conn))

se.obj <- ValidateClusters(se.obj, pc.use = 1:pcs.use, top.genes = 15, min.connectivity = min.conn, 
                           acc.cutoff = acc.cutoff, verbose = T)
se.obj <- BuildClusterTree(se.obj, do.reorder = T, reorder.numeric = T)

se.obj <- RunTSNE(se.obj, dims.use = 1:pcs.use)
TSNEPlot(se.obj)
if (regress.batch) {
  TSNEPlot(se.obj, group.by = "batch")
}

print("Done with clustering")
save.image(output.file)
saveRDS(se.obj, file = gsub(".RData", ".Robj", output.file))
