library(methods)
library(perturbLM)
library(swne)
library(glmnet)

## Command line arguments
options <- commandArgs(trailingOnly = T)

matrix.file <- options[[1]] ## Input counts matrix
genotypes.file <- options[[2]] ## Genotypes dictionary
output.file <- options[[3]] ## Output RData file
regress.batch <- as.logical(options[[4]]) ## Logical T/F whether or not there are multiple batches

## Global parameters
ctrl <- "mCherry"
alpha <- 0.5 ## Ratio of Lasso to Ridge penalty
nfolds <- 4 ## Folds for N-fold cross validations
nlambda <- 10 ## Number of shrinkage parameters to cross validate
lambda.min.ratio <- 1e-8 ## Minimum shrinkage parameter
use.quantiles <- T ## Use quantiles when calculating coefficient p-values
output.name <- 'cf' ## Name of the output variable
min.regress.cells <- 20 ## Minimum number of cells to test for marker genes
max.guides <- 2 ## Max number of guides detected per cell
n.rand <- 480 ## Number of permutations
n.cores <- 24 ## Number of cores
min.cells.frac <- 0.01 ## Minimum fraction of cells a gene must be detected in
min.genes.exp <- 150 ## Minimum number of genes a cell must express
trim <- 0.001 ## Fraction of cells to trim
n.bins <- 20 ## Number of bins to use for calculating p-values
drop.cells <- T

for (opt in options) {
  print(opt)
}

## Load the dataset, filter and trim counts matrix
counts <- ReadData(matrix.file)

## Load genotypes
genotypes.list <- ReadGenotypes(genotypes.file)
counts <- counts[,colnames(counts) %in% unlist(genotypes.list, F, F)]

# Filter for genes expressed in a minimum number of cells
counts <- FilterData(counts, min.cells.frac, trim, min.genes.exp)

## Only keep cells with a called genotype
genotypes.list <- lapply(genotypes.list, function(cells) cells[cells %in% colnames(counts)])
genotypes.list <- genotypes.list[sapply(genotypes.list, length) > 0]
print(dim(counts))

## Obtain technical covariates
print(paste("Regress out batch effects:", regress.batch))
nUMI <- scale(Matrix::colSums(counts))
if (regress.batch) {
  cell.batch <- factor(sapply(colnames(counts), function(x) strsplit(x, split = "\\.")[[1]][[2]]))
  batch.matrix <- model.matrix(~. + 0, data = data.frame(cell.batch))
  covs.df <- cbind(batch.matrix, nUMI)
  colnames(covs.df) <- c(levels(cell.batch), "nUMI")
} else {
  cell.batch <- NULL
  covs.df <- data.frame(nUMI)
}

## Build design matrix
guides.matrix <- DesignMatrixGenotypes(genotypes.list, max.guides, min.cells = 1, drop.cells = drop.cells)
guides.matrix <- CleanDesignCtrl(guides.matrix, ctrl = ctrl)

n_cells <- Matrix::colSums(guides.matrix)
guides.keep <- names(n_cells[n_cells >= min.regress.cells])
guides.matrix <- guides.matrix[,guides.keep]; n_cells <- n_cells[guides.keep];
guides.matrix <- guides.matrix[Matrix::rowSums(guides.matrix) > 0,]
print(sort(n_cells, decreasing = T))

covs.df <- data.frame(covs.df[rownames(guides.matrix),])

## Log-normalize counts matrix
log.counts <- Matrix::t(ScaleCounts(counts, batch = cell.batch, method = "log", adj.var = F))
log.counts <- log.counts[rownames(guides.matrix),]
print(dim(log.counts)); rm(counts);

## Run model fitting
cv.fit <- cv.glmnet(Matrix(cbind(as.matrix(covs.df), guides.matrix)), as.matrix(log.counts), family = "mgaussian",
                    alpha = alpha, nlambda = nlambda, nfolds = nfolds, standardize = F,
                    lambda.min.ratio = lambda.min.ratio)
lambda.use <- cv.fit$lambda.min; print(lambda.use);
save.image(output.file);

pt <- proc.time()
coefs.df <- CalcGlmnetPvals(guides.matrix, covs.df, log.counts, alpha, lambda.use, "mgaussian", 
                            ctrl = ctrl, n.rand = n.rand, n.cores = n.cores, n.bins = n.bins, 
                            use.quantiles = use.quantiles, output.name = output.name)
print(proc.time() - pt)

print("Done with regression analysis")
save.image(output.file)
write.table(coefs.df, file = gsub(".RData", ".pvals.tsv", output.file), sep = '\t')
