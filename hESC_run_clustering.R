library(methods)
library(pagoda2)
library(swne)
library(perturbLM)

options <- commandArgs(trailingOnly = T)

## Command line parameters
matrix.file <- options[[1]] ## Input counts matrix (either in tab separated format or 10X genomics sparse format)
genotypes.file <- options[[2]] ## Genotype dictionary as a ".csv" file
output.p2.file <- options[[3]] ## Output Robj file
regress.batch <- as.logical(options[[4]]) ## Whether to run batch regression

min.cells.frac <- 0.01 ## Fraction of cells expressing gene in order to include
min.genes.exp <- 200 ## Minimum number of expressed genes for a cell
trim <- 0.0025 ## Trim fraction
n.cores <- 8 ## Number of cores

# Load the dataset, filter and trim counts matrix
counts <- ReadData(matrix.file)

# Load genotypes
genotypes.list <- ReadGenotypes(genotypes.file)
genotypes.list <- lapply(genotypes.list, function(cells) cells[cells %in% colnames(counts)])
genotypes.list <- genotypes.list[sapply(genotypes.list, length) > 0]

## Filter counts
counts <- counts[,colnames(counts) %in% unlist(genotypes.list, F, F)]
counts <- FilterData(counts, min.cells.frac, trim, min.genes.exp)
print(dim(counts))

# Check for batch effects
if (regress.batch) {
  batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
} else {
  batch <- NULL
}

## Create pagoda2 object
rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts, batch = batch, log.scale = T, n.cores = n.cores)

## Adjust variance, calculate PCA
r$adjustVariance(plot = T, gam.k = 10)
r$calculatePcaReduction(nPcs = 50, n.odgenes = 3e3)

## kNN clustering
r$makeKnnGraph(k = 50, type='PCA', center = T, distance = 'cosine')
r$getKnnClusters(method = multilevel.community, type = 'PCA', name = "multilevel")

## tSNE
r$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 50, verbose = F)

## Differential expression
r$getDifferentialGenes(type = 'PCA', verbose = T, clusterType = 'multilevel')

## Save output
saveRDS(r, file = output.p2.file)
