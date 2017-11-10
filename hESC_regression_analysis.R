library(methods)
library(snow)
library(Seurat)
source('hESC_functions.R')

# Command line arguments
options <- commandArgs(trailingOnly = T)

matrix.file <- options[[1]] ## Input counts matrix
genotypes.file <- options[[2]] ## Genotypes dictionary
output.file <- options[[3]] ## Output RData file
regress.batch <- as.logical(options[[4]]) ## Logical T/F whether or not there are multiple batches

# Global parameters
alpha <- 0.5 ## Ratio of Lasso to Ridge penalty
nfolds <- 5 ## Folds for N-fold cross validations
nlambda <- 20 ## Number of shrinkage parameters to cross validate
lambda.min.ratio <- 1e-10 ## Minimum shrinkage parameter
use.quantiles <- T ## Use quantiles when calculating coefficient p-values
output.name <- 'cf' ## Name of the output variable
min.regress.cells <- 25 ## Minimum number of cells to test for marker genes
max.guides <- 2 ## Max number of guides detected per cell
n.rand <- 200 ## Number of permutations
n.cores <- 8 ## Number of cores
min.cells.frac <- 0.025 ## Minimum fraction of cells a gene must be detected in
min.genes.exp <- 200 ## Minimum number of genes a cell must express
trim <- 0.005 ## Fraction of cells to trim
n.bins <- 20 ## Number of bins to use for calculating p-values
drop.cells <- F

for (opt in options) {
  print(opt)
}

# Load the dataset, filter and trim counts matrix
counts <- read.data(matrix.file)

# Load genotypes
genotypes.list <- read.genotypes(genotypes.file)
counts <- counts[,colnames(counts) %in% unlist(genotypes.list, F, F)]

# Filter for genes expressed in a minimum number of cells
counts <- filter.data(counts, min.cells.frac = min.cells.frac, trim = trim, min.genes = min.genes.exp, min.expr = 0)

# Only keep cells with a called genotype
genotypes.list <- lapply(genotypes.list, function(cells) cells[cells %in% colnames(counts)])
genotypes.list <- genotypes.list[sapply(genotypes.list, length) > 0]
print(dim(counts))

# Obtain technical covariates
print(paste("Regress out batch effects:", regress.batch))
nUMI <- scale(colSums(counts))
if (regress.batch) {
  cell.batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
  batch.matrix <- model.matrix(~. + 0, data = data.frame(cell.batch))
  covs.df <- cbind(batch.matrix, nUMI)
  colnames(covs.df) <- c(levels(cell.batch), "nUMI")
} else {
  covs.df <- data.frame(nUMI)
}

# Build design matrix
guides.matrix <- Matrix(design.matrix.genotypes(genotypes.list, max.guides, min.regress.cells, drop.cells = drop.cells))
n_cells <- Matrix::colSums(guides.matrix)
covs.df <- covs.df[rownames(guides.matrix),]
print(n_cells)

# Log-normalize counts matrix
log.counts <- t(log(counts + 1))
log.counts <- log.counts[rownames(guides.matrix),]
rm(counts)

# Run model fitting
cv.fit <- cv.glmnet(cbind(covs.df, guides.matrix), y = as.matrix(log.counts), family = "mgaussian", alpha = alpha,
                    nlambda = nlambda, nfolds = nfolds, standardize = F, lambda.min.ratio = lambda.min.ratio )
best.lambda <- cv.fit$lambda.min
plot(cv.fit)
print(best.lambda)

coefs.df <- calc.glmnet.pvals(guides.matrix, covs.df, log.counts, alpha, best.lambda, "mgaussian", n.rand, n.cores, binned.pvals = T,
                              n.bins = n.bins, use.quantiles = use.quantiles, output.name = output.name)

print("Done with genotype enrichment")
save.image(output.file)
write.table(coefs.df, file = gsub(".RData", ".pvals.tsv", output.file), sep = '\t')

