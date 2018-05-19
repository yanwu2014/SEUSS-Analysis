library(perturbLM)

options <- commandArgs(trailingOnly = T)
ctrl <- "mCherry"
seed <- 312513
frac <- 0.25

## Load genotypes
genotypes.file <- options[[1]]
genotypes.list <- ReadGenotypes(genotypes.file)

## Create internal control
set.seed(seed)
ctrl.cells <- genotypes.list[[ctrl]]
int.ctrl.cells <- sample(ctrl.cells, round(frac*length(ctrl.cells)))
genotypes.list[[paste(ctrl, "int-ctrl", sep = "-")]] <- int.ctrl.cells
genotypes.list[[ctrl]] <- ctrl.cells[!ctrl.cells %in% int.ctrl.cells]

WriteGenotypes(genotypes.list, genotypes.file)