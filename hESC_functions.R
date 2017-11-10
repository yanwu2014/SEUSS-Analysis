## Helper functions ##
library(matrixStats)
library(compiler)
library(reshape2)
library(Matrix)
library(MASS)
library(abind)
library(glmnet)
library(lsr)
library(reshape2)
library(data.table)


# Basic function to convert vector of mouse genes to human orthologs
convertMouseGeneList <- function(x) {
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,2])
  return(humanx)
}


## Trim input matrices. From PAGODA2
winsorize.matrix <- function(mat, trim) {
  require(Rcpp)
  sourceCpp("hESC_cpp_functions.cpp")
  
  if(trim  >  0.5) { trim <- trim/ncol(mat)  }
  wm <- winsorizeMatrix(mat, trim)
  rownames(wm) <- rownames(mat)
  colnames(wm) <- colnames(mat)
  return(wm)
}

winsorize.matrix <- cmpfun(winsorize.matrix)


## Build Shared Nearest Neighbors matrix. From Seurat
BuildSNN <- function(data.use, k.param = 10, k.scale = 10, prune.SNN = 1/15) {
  require(FNN)
  n.cells <- nrow(data.use)
  
  my.knn <- get.knn(data = data.use, k = min(k.scale * k.param, n.cells - 1))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param - 1)])
  nn.large <- my.knn$nn.index
  
  w <- CalcSNNSparse(cell.names = rownames(data.use), k.param = k.param, nn.large = nn.large,
                     nn.ranked = nn.ranked, prune.SNN = prune.SNN)
  
  return(w)
}

BuildSNN <- cmpfun(BuildSNN)


## Helper function for building SNN from Seurat
CalcSNNSparse <- function(cell.names, k.param, nn.large, nn.ranked, prune.SNN, print.output = T) {
  n.cells <- length(cell.names)
  counter <- 1
  idx1 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  idx2 <- vector(mode = "integer", length = n.cells ^ 2 / k.param)
  edge.weight <- vector(mode = "double", length = n.cells ^ 2 / k.param)
  id <- 1
  # fill out the adjacency matrix w with edge weights only between your target
  # cell and its k.scale*k.param-nearest neighbors
  # speed things up (don't have to calculate all pairwise distances)
  # define the edge weights with Jaccard distance
  if (print.output) {
    print("Constructing SNN")
    pb <- txtProgressBar(min = 0, max = n.cells, style = 3)
  }
  for (i in 1:n.cells) {
    for (j in 1:ncol(x = nn.large)) {
      s <- intersect(x = nn.ranked[i, ], y = nn.ranked[nn.large[i, j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i, j], ])
      e <- length(x = s) / length(x = u)
      if (e > prune.SNN) {
        idx1[id] <- i
        idx2[id] <- nn.large[i, j]
        edge.weight[id] <- e
        id <- id + 1
      }
    }
    if (print.output) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (print.output) {
    close(con = pb)
  }
  idx1 <- idx1[! is.na(x = idx1) & idx1 != 0]
  idx2 <- idx2[! is.na(x = idx2) & idx2 != 0]
  edge.weight <- edge.weight[! is.na(x = edge.weight) & edge.weight != 0]
  w <- sparseMatrix(
    i = idx1,
    j = idx2,
    x = edge.weight,
    dims = c(n.cells, n.cells)
  )
  diag(x = w) <- 1
  rownames(x = w) <- cell.names
  colnames(x = w) <- cell.names
  return(w)
}

CalcSNNSparse <- cmpfun(CalcSNNSparse)


## Calculate connectivity between clusters given an SNN matrix. From Seurat
CalcConnectivity <- function(object) {
  SNN <- object@snn
  cluster.names <- unique(x = object@ident)
  num.clusters <- length(x = cluster.names)
  connectivity <- matrix(data = 0, nrow = num.clusters, ncol = num.clusters)
  rownames(x = connectivity) <- cluster.names
  colnames(x = connectivity) <- cluster.names
  n <- 1
  for (i in cluster.names) {
    for (j in cluster.names[-(1:n)]) {
      subSNN <- SNN[
        match(x = WhichCells(object = object, ident = i), colnames(x = SNN)),
        match(x = WhichCells(object = object, ident = j), rownames(x = SNN))
        ]
      if (is.object(x = subSNN)) {
        connectivity[i, j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[i, j] <- mean(x = subSNN)
      }
    }
    n <- n + 1
  }
  return(connectivity)
}


#### Functions for handling genotype/phenotype dictionaries #####

## Convert genotypes dictionary from list to named vector. All cells with more than one genotype are removed
flatten.genotype.list <- function(genotypes.list) {
  cell.names <- unlist(genotypes.list, F, F)
  if (length(cell.names) > length(unique(cell.names))) {
    print("Warning: removing all cells with duals")
    genotypes.list <- single.genotypes(genotypes.list)
    cell.names <- unlist(genotypes.list, F, F)
  }
  
  genotypes <- rep(NA, length(cell.names))
  names(genotypes) <- cell.names
  for (g in names(genotypes.list)) {
    genotypes[genotypes.list[[g]]] <- g
  }
  
  if (any(is.na(genotypes))) { stop("Some unassigned cells"); }
  return(genotypes)
}


single.genotypes <- function(genotypes.list, min.cells = 30) {
  genotypes.matrix <- design.matrix.genotypes(genotypes.list, max.guides = 1, min.cells = min.cells)
  genotypes.list <- lapply(colnames(genotypes.matrix), function(g) {
    rownames(genotypes.matrix)[which(genotypes.matrix[,g] == 1)]
  })
  names(genotypes.list) <- colnames(genotypes.matrix)
  return(genotypes.list)
}

single.genotypes <- cmpfun(single.genotypes)


## Convert genotypes from named vector to list
unflatten.cell.genotypes <- function(cell.genotypes, min.cells = 1) {
  genotypes.list <- c()
  genotypes <- unique(cell.genotypes)
  for (g in genotypes) {
    g.cells <- names(cell.genotypes[cell.genotypes == g])
    if (length(g.cells) > min.cells) {
      genotypes.list[[g]] <- g.cells
    }
  }
  return(genotypes.list)
}

unflatten.cell.genotypes <- cmpfun(unflatten.cell.genotypes)


## Read in genotypes dictionary from csv file. Returns an R list.
read.genotypes <- function(pheno.dict, sep.char = ",") {
  geno.data <- read.table(pheno.dict, sep = sep.char, header = F, stringsAsFactors = F)
  genotypes <- geno.data[[1]]
  full.genotypes.list <- lapply(rownames(geno.data), function(i) sapply(strsplit(geno.data[i,2], split = ',')[[1]], trimws))
  full.genotypes.list <- lapply(full.genotypes.list, function(x) make.names(x))
  names(full.genotypes.list) <- genotypes
  return(full.genotypes.list)
}

read.genotypes <- cmpfun(read.genotypes)


## Write R list to csv file
write.genotypes <- function(genotypes.list, out.file) {
  geno.data <- sapply(genotypes.list, function(x) paste('\"', paste(x, collapse = ", "), '\"', sep = ""))
  geno.data <- sapply(names(geno.data), function(x) paste(x, geno.data[[x]], sep = ","))
  fileConn = file(out.file)
  writeLines(geno.data, fileConn)
  close(fileConn)
}

write.genotypes <- cmpfun(write.genotypes)


## Convert genotypes list to design matrix
design.matrix.genotypes <- function(full.genotypes.list, max.guides = 2, min.cells = 5, drop.cells = T) {
  require(hash)
  cell.names <- unique(unlist(full.genotypes.list, F, F))
  single.genotypes <- names(full.genotypes.list)
  
  single.mat <- Matrix(0, length(cell.names), length(single.genotypes), dimnames = list(cell.names, single.genotypes))
  for (i in 1:ncol(single.mat)) {
    single.mat[full.genotypes.list[[i]], i] <- 1
  }
  if (drop.cells) single.mat <- single.mat[Matrix::rowSums(single.mat) <= max.guides, ]
  n_cells_single <- Matrix::colSums(single.mat)
  single.mat <- single.mat[ ,names(n_cells_single[n_cells_single >= min.cells])]
  single.mat <- single.mat[Matrix::rowSums(single.mat) > 0, ]
  
  combo.rows <- single.mat[Matrix::rowSums(single.mat) > 1,]
  combo.genotypes <- unique(apply(combo.rows, 1, function(x) paste(names(x[x == 1]), collapse = ":", sep = "")))
  if (max.guides > 1 && length(combo.genotypes) > 0) {
    combo.genotypes.list <- hash(keys = combo.genotypes)
    for (i in 1:nrow(combo.rows)) {
      mat.row <- combo.rows[i,]
      genotype <- paste(names(mat.row[mat.row == 1]), collapse = ":", sep = "")
      cell <- rownames(combo.rows)[[i]]
      combo.genotypes.list[[genotype]] <- c(combo.genotypes.list[[genotype]], cell)
    }
    combo.genotypes.list <- as.list(combo.genotypes.list)
    combo.genotypes.list[['keys']] <- NULL
    combo.genotypes.list <- combo.genotypes.list[combo.genotypes]
    
    combo.mat <- Matrix(0, nrow(single.mat), length(combo.genotypes), dimnames = list(rownames(single.mat), combo.genotypes))
    for (i in 1:ncol(combo.mat)) {
      combo.mat[combo.genotypes.list[[i]],i] <- 1
    }
    
    # Filter out combos that have too few cells    
    n_cells_combo <- Matrix::colSums(combo.mat)
    combos_keep <- names(n_cells_combo[n_cells_combo >= min.cells])
    
    # Remove filtered combo cells from design matrix
    if (length(combos_keep) > 0) {
      combo.mat <- combo.mat[ ,combos_keep]
      design.mat <- cbind(single.mat, combo.mat)
      filtered.cells <- (Matrix::rowSums(single.mat) > 1) & (Matrix::rowSums(combo.mat) == 0)
      
    } else {
      design.mat <- single.mat
      filtered.cells <- (Matrix::rowSums(single.mat) > 1)
      
    }
    if (drop.cells) design.mat <- design.mat[!filtered.cells,]
  } else {
    design.mat <- single.mat
  }
  return(design.mat)
}

design.matrix.genotypes <- cmpfun(design.matrix.genotypes)


## Pad design matrix with zeros for cells without a called genotype
pad.design.matrix <- function(X, row.names) {
  Y <- matrix(0, length(row.names), ncol(X), dimnames = list(row.names, colnames(X)))
  Y[rownames(X),] <- X
  return(Y)
}

pad.design.matrix <- cmpfun(pad.design.matrix)


## Given a list of cell clusters and a list of cell genotypes,
## Returns a genotypes x clusters matrix of cell counts
genotype.cluster.counts <- function(genotypes.list, clusters.list) {
  
  df <- matrix(0, length(genotypes.list), length(clusters.list))
  rownames(df) <- names(genotypes.list)
  colnames(df) <- names(clusters.list)
  
  for(i in 1:length(genotypes.list)) {
    genotype.cells <- genotypes.list[[i]]
    for(j in 1:length(clusters.list)) {
      cluster.cells <- clusters.list[[j]]
      df[i,j] <- length(intersect(genotype.cells, cluster.cells))
    }
  }
  return(df)
}


#### Handling counts matrices ####

## Read data from either a tab separated values file or a 10X genomics sparse matrix directory
read.data <- function(matrix.dir) {
  if (dir.exists(matrix.dir)) {
    counts <- as.matrix(Seurat::Read10X(matrix.dir))
  } else {
    counts <- read.table(matrix.dir, sep = "\t", header = T, row.names = 1)
  }
  colnames(counts) <- make.names(colnames(counts))
  return(counts)
}

read.data <- cmpfun(read.data)


# Filter data by minimum cells, minimum genes
filter.data <- function(x, min.cells.frac, trim, min.genes = 500, min.expr = 0, max.transcripts = 70000) {
  min.cells <- round(ncol(x)*min.cells.frac)
  x <- x[ , Matrix::colSums(x) < max.transcripts]
  x <- x[ , Matrix::colSums(x > min.expr) > min.genes]
  x <- x[Matrix::rowSums(x > min.expr) > min.cells, ]
  if (trim > 0) {
    x <- t(winsorize.matrix(t(x), trim = trim))
  }
  return(x)
}

filter.data <- cmpfun(filter.data)


#### Functions for handling genesets ####

## Read genesets from gmt file to list format
load.genesets <- function(geneset.name) {
  geneset.to.use <- as.list(readLines(geneset.name))
  geneset.to.use <- lapply(geneset.to.use, function (v) strsplit(v, '\t')[[1]])
  geneset.names <- unlist(lapply(geneset.to.use, function(x) x[[1]]))
  geneset.to.use <- lapply(geneset.to.use, function(v) v[3:length(v)])
  names(geneset.to.use) <- geneset.names
  names(geneset.to.use) <- sapply(names(geneset.to.use), function(x) gsub(" ", "_", x))
  return(geneset.to.use)
}

## Filter genesets by minimum, maximum number of genes
filter.genesets <- function(geneset.to.use, gene.names, min.size = 5, max.size = 500) {
  geneset.to.use <- lapply(geneset.to.use, function (x) return(x[x %in% gene.names]))
  geneset.to.use <- clean.genesets(geneset.to.use, min.size = min.size, max.size = max.size)
  return(geneset.to.use)
}


## Helper function for filtering genesets
clean.genesets <- function(go.env, min.size = 5, max.size = 500, annot = FALSE) {
  go.env <- as.list(go.env)
  size <- unlist(lapply(go.env, length))
  go.env <- go.env[size > min.size & size < max.size]
  return(go.env)
}


## Write genesets from list to gmt file
write.genesets <- function(genesets, file.name) {
  genesets <- lapply(names(genesets), function(name) {x <- genesets[[name]]; x <- c(name,name,x); return(x);})
  n.cols <- 1.5*max(unlist(lapply(genesets, length)))
  empty <- lapply(genesets, write, file.name, append=TRUE, ncolumns = n.cols, sep = '\t')
}


calc.gene_modules.effect <- function(coefs.matrix, gene_clusters_list, gene_module_mapping = NULL, min.coef = 0.025) {
  coefs.matrix <- coefs.matrix[apply(coefs.matrix, 1, function(x) any(abs(x) > min.coef)),]
  
  gene_clusters_list <- load.genesets("up-tf_functional_gene_modules.gmt")
  if (!is.null(gene_module_mapping)) {
    gene_modules <- gene_module_mapping$Description
    names(gene_modules) <- gene_module_mapping$Module
    names(gene_clusters_list) <- sapply(names(gene_clusters_list), function(x) gene_modules[[x]])
  }
  gene_clusters <- flatten.genotype.list(gene_clusters_list)
  
  coefs.matrix <- coefs.matrix[rownames(coefs.matrix) %in% names(gene_clusters),]
  gene_clusters <- gene_clusters[rownames(coefs.matrix)]
  
  tfs.modules.matrix <- apply(coefs.matrix, 2, function(x) {
    tapply(x, gene_clusters, mean)
  })
  tfs.modules.matrix <- tfs.modules.matrix[,order(colnames(tfs.modules.matrix))]
  tfs.modules.matrix
}


#### Functions for flattening and unflattening dataframes #####

## Convert dataframe to matrix, specifying all column names
unflatten.dataframe <- function(df, output.name, row.col = 'Gene', col.col = 'Group', 
                                direction = F) {
  if (direction) {
    df$Gene <- paste(df$Gene, df$Direction, sep = "_")
  }
  df <- df[c(row.col, col.col, output.name)]
  colnames(df) <- c('row', 'column', output.name)
  mat.out <- acast(df, row~column, value.var = output.name)
  return(mat.out)
}

unflatten.dataframe <- cmpfun(unflatten.dataframe)


#### Visualization functions ####

## Visualize two dimensions
dim_plot <- function(dim.x, dim.y, color.plot = NULL, x.lab = "tsne1", y.lab = "tsne2", main.title = NULL,
                      pt.size = 1.0, font.size = 12, alpha = 1.0, group.label = T, label.size = 4, frac.trim = 0.025,
                      show.legend = T, label.points = F) {
  require(ggplot2)
  require(ggrepel)
  gg.df <- data.frame(x = dim.x, y = dim.y, color.plot = color.plot, pts.label = names(dim.x))
  ggobj <- ggplot(gg.df, aes(x, y)) + geom_point(size = pt.size, alpha = alpha, aes(colour = color.plot)) +
    theme_void() + theme(text = element_text(size = font.size)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title)
  
  if (group.label && is.factor(color.plot)) {
    group.pts.x <- tapply(gg.df$x, gg.df$color.plot, function(v) {
      v <- DescTools::Trim(sort(v), trim = frac.trim)  
      median(v)
    })
    group.pts.y <- tapply(gg.df$y, gg.df$color.plot, function(v) {
      v <- DescTools::Trim(sort(v), trim = frac.trim)
      median(v)
    })
    group.pts <- data.frame(x = group.pts.x, y = group.pts.y)
    group.pts$ident <- levels(gg.df$color.plot)
    
    ggobj <- ggobj + geom_point(data = group.pts, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
      geom_text_repel(data = group.pts, mapping = aes(label = ident), size = label.size)
  }
  
  if (!show.legend) {
    ggobj <- ggobj + theme(legend.position = "none")
  }
  
  if (label.points) {
    ggobj <- ggobj + geom_text_repel(mapping = aes(label = pts.label), size = label.size)
  }
  
  ggobj
}


## Plot matrix as heatmap. Options for clustering rows/columns
ggheat <- function(m, rescaling = 'none', clustering = 'none', labCol = T, labRow = T, border = F, 
                   heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato'), legend.title = NULL, x.lab.size = 8,
                   y.lab.size = 8) {
  require(reshape)
  require(ggplot2)
  
  ## you can either scale by row or column not both! 
  ## if you wish to scale by both or use a differen scale method then simply supply a scale
  ## function instead NB scale is a base funct
  
  if(is.function(rescaling)) { 
    m = rescaling(m)
  } else 
  {
    if(rescaling == 'column') 
      m = scale(m, center=T)
    if(rescaling == 'row') 
      m = t(scale(t(m),center=T))
  }
  
  ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function 
  
  if(is.function(clustering)) {
    m = clustering(m)
  } else {
    if(clustering=='row')
      m = m[hclust(dist(m))$order, ]
    if(clustering=='column')  
      m = m[ ,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m = m[hclust(dist(m))$order, hclust(dist(t(m)))$order]
  }
  ## this is just reshaping into a ggplot format matrix and making a ggplot layer
  
  rows = dim(m)[1]
  cols = dim(m)[2]
  melt.m = cbind(rowInd=rep(1:rows, times = cols), colInd = rep(1:cols, each = rows), melt(m))
  g = ggplot(data = melt.m)
  
  ## add the heat tiles with or without a white border for clarity
  if(border == TRUE)
    g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd, fill = value), colour = 'white')
  if(border == FALSE)
    g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd, fill = value))
  
  ## add axis labels either supplied or from the colnames rownames of the matrix
  if(labCol == T)
    g2 = g2 + scale_x_continuous(breaks = (1:cols) - 0.5, labels = colnames(m), expand = c(0.005,0))
  if(labCol == F) 
    g2 = g2 + scale_x_continuous(breaks = (1:cols) - 0.5, labels = rep('', cols))
  if(labRow == T)
    g2 = g2 + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rownames(m), expand = c(0.005,0))	
  if(labRow == F) 
    g2 = g2 + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rep('', rows))
  
  ## get rid of grey panel background and gridlines
  g2 = g2 + theme(panel.grid.minor = element_line(colour = NA), panel.grid.major = element_line(colour=NA),
                  panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_text(angle = 90, hjust = 1, size = x.lab.size),
                  axis.ticks = element_blank(), axis.text.y = element_text(size = y.lab.size))
  
  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  return(g2 + scale_fill_gradient2(low = heatscale[1], mid = heatscale[2], high = heatscale[3], guide = guide_colorbar(title = legend.title)))
}


## Visualize correlation between two vectors
ggcorrelation <- function(x, y, x.lab, y.lab, title, genes.use = NULL, use.label = F, 
                          points.label = NULL, font.size = 14, label.font.size = 4, pt.size = 1,
                          show.corr = F, box = T, alpha = 1) {
  require(ggplot2)
  require(ggrepel)
  
  stopifnot(names(x) == names(y))
  if (is.null(genes.use)) {
    genes.use <- names(x)
  }
  corr.use <- paste("R = ", round(cor(x[genes.use], y[genes.use]), 2))
  
  gg.df <- data.frame(x[genes.use],y[genes.use])
  gg.df$magnitude <- abs(x[genes.use]) + abs(y[genes.use])
  gg.df$name <- names(x[genes.use])
  
  if (use.label) {
    gg.df$label <- points.label
    gg.df$name[!points.label] <- ""
    colnames(gg.df) <- c("x", "y", 'magnitude', 'name', 'label')
  } else {
    colnames(gg.df) <- c("x", "y", 'magnitude', 'name')
  }
  
  if (show.corr) {
    main.title <- paste(title, corr.use, sep = ": ")
  } else {
    main.title <- title
  }
  
  min.pt <- min(c(min(x[genes.use]), min(y[genes.use])))
  max.pt <- max(c(max(x[genes.use]), max(y[genes.use])))
  
  ggobj <- ggplot(gg.df, aes(x, y)) + geom_point(aes(colour = magnitude), size = pt.size, alpha = alpha) +
    scale_colour_gradient(low = 'grey', high = 'blue') + theme_classic() +
    theme(legend.position="none", text = element_text(size = font.size)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title)
  if (box) {
    ggobj <- ggobj + xlim(c(min.pt, max.pt)) + ylim(c(min.pt, max.pt))
  }
  if (use.label) {
    ggobj <- ggobj + geom_text_repel(aes(x, y, label = name), size = label.font.size)
  }
  return(ggobj)
}

ggcorrelation <- cmpfun(ggcorrelation)


gghexbin <- function(v, w, n.bins = 100) {
  ii <- intersect(names(v), names(w))
  r <- round(cor(v[ii], w[ii]), 3)
  df <- data.frame(x = v[ii], y = w[ii])
  ggobj <- ggplot(df) + geom_hex(aes(x, y), bins = n.bins, alpha = 1) + 
    theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") + 
    scale_fill_gradient(low = "skyblue", high = "tomato") + ggtitle(paste("R =", r))
  ggobj
}


### Enrichment functions ####

## Use Fisher's exact test to calculate enrichment for count data
fisher.enrich <- function(marker.genes, genesets, n.background, correct = F) {
  results.df <- data.frame(genesets = names(genesets), p.val = numeric(length(genesets)),
                           genes_in_set = numeric(length(genesets)), geneset_size = numeric(length(genesets)),
                           stringsAsFactors = F)
  
  for (i in 1:length(genesets)) {
    set.genes <- genesets[[i]]
    set.size <- length(set.genes)
    in.set <- length(intersect(set.genes, marker.genes)) - 1
    out.set <- length(marker.genes) - in.set
    
    if (in.set < 0) { in.set <- 0 }
    cont.table <- matrix(c(in.set, set.size, out.set, n.background - set.size), 
                         nrow = 2, byrow = T)
    
    results.df[i, "p.val"] <- fisher.test(cont.table, alternative = "greater")$p.value
    results.df[i, "geneset_size"] <- set.size
    results.df[i, "genes_in_set"] <- in.set
  }
  
  if (correct) {
    results.df$FDR <- p.adjust(results.df$p.val, method = "BH")
  }
  
  results.df <- results.df[order(results.df$p.val),]
  return(results.df)
}

fisher.enrich <- cmpfun(fisher.enrich)


## Count overenrichment for a list of gene markers
multiple.fisher.enrich <- function(markers.list, genesets) {
  all.genes <- unique(unlist(markers.list, F, F))
  genesets <- filter.genesets(genesets, all.genes, min.size = 1, max.size = 10000)
  
  enrich.list <- lapply(names(markers.list), function(g) {
    res.df <- fisher.enrich(markers.list[[g]], genesets, length(all.genes), correct = T)
    res.df$Group <- g; rownames(res.df) <- NULL;
    res.df
  })
  enrich.df <- do.call("rbind", enrich.list)
  # enrich.df$FDR <- p.adjust(enrich.df$p.val, method = "BY")
  enrich.df[order(enrich.df$p.val),]
}

multiple.fisher.enrich <- cmpfun(multiple.fisher.enrich)


## GSEA over-enrichment test using liger
multiple.gsea.enrich <- function(scores.list, genesets, n.rand = 1000, n.cores = 1, power = 1) {
  require(liger)
  
  enrich.list <- lapply(names(scores.list), function(g) {
    res <- bulk.gsea(scores.list[[g]], genesets, mc.cores = n.cores, n.rand = n.rand, power = power, rank = F,
                     skip.qval.estimation = T)
    res$Group <- g; res$genesets <- rownames(res); res$q.val <- NULL;
    res
  })
  enrich.df <- do.call("rbind", enrich.list); rownames(enrich.df) <- NULL;
  enrich.df$FDR <- p.adjust(enrich.df$p.val, method = "BY")
  enrich.df[order(enrich.df$p.val),]
}

multiple.gsea.enrich <- cmpfun(multiple.gsea.enrich)


## Calculate chi-squared test of difference between each row of a counts matrix and the specified control row
genotype.cluster.chisq <- function(df, ntc) {
  apply(df, 1, function(x) {
    chisq.test(rbind(x, df[ntc,]))$p.value
  })
}


## Use Fisher's exact test to calculate enrichment for each genotype in each cluster
genotype.cluster.pvals <- function(df.counts) {
  cluster.counts <- colSums(df.counts)
  genotype.counts <- rowSums(df.counts)
  n.cells <- sum(cluster.counts)
  
  p.vals <- matrix(1, nrow = nrow(df.counts), ncol = ncol(df.counts),
                   dimnames = list(rownames(df.counts), colnames(df.counts)))
  
  for(i in 1:nrow(df.counts)) {
    for(j in 1:ncol(df.counts)) {
      x <- df.counts[i,j] - 1
      if (x < 0) { x <- 0; }
      cont.table <- matrix(c(x, cluster.counts[[j]], genotype.counts[[i]] - x, n.cells), nrow = 2, byrow = T)
      p.vals[i,j] <- fisher.test(cont.table, alternative = "greater")$p.value
    }
  }
  
  return(p.vals)
}


#### Linear models ####

## Extract coefficient matrix from a multigaussian glment object
get.coef.matrix <- function(mfit, best.lambda) {
  cfs <- coef(mfit, s = best.lambda)
  cfs <- lapply(cfs, function(x) {y <- as.numeric(x); names(y) <- rownames(x); return(y)})
  cfs <- do.call('rbind', cfs)
  colnames(cfs) <- gsub('genotype', '', colnames(cfs))
  return(cfs)
}

get.coef.matrix <- cmpfun(get.coef.matrix)


## Calculate regression coefficients and p-values via permutation testing
calc.glmnet.pvals <- function(design.matrix, metadata, y, alpha, lambda, family, n.rand, n.cores, binned.pvals = F,
                              n.bins = 10, use.quantiles = T, output.name = "cf") {
  mfit <- glmnet(Matrix(cbind(design.matrix, metadata)), y = y, family = family, alpha = alpha, lambda = lambda,
                 standardize = F)
  cfs <- get.coef.matrix(mfit, lambda)
  
  cols.keep <- colnames(cfs)[2:ncol(cfs)]
  cols.keep <- cols.keep[!cols.keep %in% colnames(metadata)]
  
  cl <- makeCluster(n.cores, type = "SOCK")
  clusterExport(cl, c("design.matrix", "y", "metadata", "alpha", "lambda", "family", "get.coef.matrix"), envir = environment())
  source.log <- parLapply(cl, 1:n.cores, function(i) library(glmnet))
  cfs.rand <- parLapply(cl, 1:n.rand, function(i) {
    design.matrix.permute <- design.matrix[sample(1:nrow(design.matrix)),]
    mfit <- glmnet(Matrix(cbind(design.matrix.permute, metadata)), y = y, family = family, alpha = alpha, lambda = lambda,
                   standardize = F)
    get.coef.matrix(mfit, lambda)
  })
  
  if (binned.pvals) {
    gene_vars <- matrixStats::colVars(y)
    gene_avgs <- colMeans(y)
    gene.covs <- data.frame(gene_avgs, gene_vars)
    
    cfs <- cfs[,cols.keep]
    cfs.rand <- lapply(cfs.rand, function(x) x[,cols.keep])
    
    covs.use <- colnames(gene.covs)
    # coefs.df <- calc.binned.pvals(cfs, cfs.rand, gene.covs, NULL, covs.use, n.bins, use.quantiles, output.name, round(n.cores/4))
    
    null.coefs.df <- lapply(cfs.rand, flatten.score.matrix, output.name = output.name, gene.covs = gene.covs, group.covs = NULL)
    null.coefs.df <- rbindlist(null.coefs.df)
    null.coefs.df <- get.multi.bins(null.coefs.df, covs.use, n.bins, use.quantiles, bin.group = T)
    null.binned.dat <- get.binned.list(null.coefs.df, output.name)
    rm(null.coefs.df)
    gc()
    
    coefs.df <- flatten.score.matrix(cfs, output.name, gene.covs, NULL)
    coefs.df <- get.multi.bins(coefs.df, covs.use, n.bins, use.quantiles, bin.group = T)
    coefs.df <- calc.emp.pvals(coefs.df, null.binned.dat, output.name = 'cf', n.cores = round(n.cores/4))
    coefs.df <- coefs.df[order(coefs.df$p_val),]
    
    coefs.df[c("gene_avgs", "gene_vars", "bin.index")] <- NULL
    return(coefs.df)
  } else {
    cfs.rand <- abind::abind(cfs.rand, along = 3)
    cfs.pvals <- calc.pvals.cfs(cfs, cfs.rand)
    
    fit.res <- list()
    fit.res$coefs <- cfs[,cols.keep]
    fit.res$pvals <- cfs.pvals[,cols.keep]
    return(fit.res)
  }
}

calc.glmnet.pvals <- cmpfun(calc.glmnet.pvals)



#### Permutation testing and P-value matrix manipulation ####

## Given a matrix of regression coefficients, and a 3D array of permuted coefficients
## calculate empirical p-values
calc.pvals.cfs <- function(cfs, cfs.rand) {
  p.vals <- matrix(1, nrow(cfs), ncol(cfs))
  colnames(p.vals) <- colnames(cfs)
  rownames(p.vals) <- rownames(cfs)
  for (i in 1:nrow(cfs)) {
    for (j in 1:ncol(cfs)) {
      v <- sort(na.omit(cfs.rand[i,j,]))
      v.pos <- v[v > 0]
      v.neg <- v[v <= 0]
      b <- cfs[i,j]
      if (b > 0) {
        p.vals[i,j] <- calcPvalGreaterCpp(v.pos, b)
      } else {
        p.vals[i,j] <- calcPvalLessCpp(v.neg, b)
      }
    }
  }
  return(p.vals)
}

calc.pvals.cfs <- cmpfun(calc.pvals.cfs)


## Convert matrix of p-values into a dataframe
flatten.pvals <- function(pvals, effect, bi = "F", correct = F) {
  pvals <- pvals[rownames(effect),]
  stopifnot(all(rownames(pvals) == rownames(effect)))
  if (bi) {
    stopifnot(ncol(pvals) == 2*ncol(effect))
  } else {
    stopifnot(ncol(pvals) == ncol(effect) && all(colnames(pvals) == colnames(effect)))
  }
  
  # Organize top hits into data frame
  num.rows <- length(pvals)
  top.hits.df <- data.frame(Group = character(num.rows), Gene = character(num.rows),
                            p_val = numeric(num.rows), log.fc = numeric(num.rows),
                            stringsAsFactors = F)
  i <- 1
  for (genotype in rownames(pvals)) {
    p.vals <- pvals[genotype,]
    ef.sizes <- effect[genotype,]
    features <- names(p.vals)
    
    j <- i + length(p.vals) - 1
    top.hits.df[i:j,1] <- genotype
    top.hits.df[i:j,2] <- features
    top.hits.df[i:j,3] <- p.vals
    
    if (bi) {
      ef.sizes <- c(ef.sizes, ef.sizes)
    }
    top.hits.df[i:j,4] <- ef.sizes
    i <- j + 1
  }
  
  top.hits.df <- top.hits.df[order(top.hits.df$p_val), ]
  if (correct) {
    top.hits.df$FDR <- qvalue::qvalue(top.hits.df$p_val)$qvalues
  }
  return(top.hits.df)
}

flatten.pvals <- cmpfun(flatten.pvals)


## Helper function that converts a list to a matrix
list.to.matrix <- function(l) {
  max_length <- max(sapply(l,length))
  vapply(l, function(x){
    c(x, rep(NA, max_length - length(x)))
  }, rep(1.0, max_length))
}

list.to.matrix <- cmpfun(list.to.matrix)


## Helper function for p-value calculation. 
flatten.score.matrix <- function(score.matrix, output.name, gene.covs = NULL, group.covs = NULL, add.direction = F) {
  scores.df <- as.data.frame.table(score.matrix, stringsAsFactors = F, responseName = output.name)
  colnames(scores.df) <- c('Gene','Group', output.name)
  num.cfs <- nrow(scores.df)
  
  if (!is.null(gene.covs)) {
    for (cov in colnames(gene.covs)) {
      scores.df[[cov]] <- numeric(num.cfs)
    }
  }
  
  if (!is.null(group.covs)) {
    for (cov in colnames(group.covs)) {
      scores.df[[cov]] <- numeric(num.cfs)
    }
  }
  
  st <- 1
  en <- nrow(score.matrix)
  for (i in 1:ncol(score.matrix)) {
    for (cov in colnames(gene.covs)) {
      scores.df[st:en, cov] <- gene.covs[[cov]]
    }
    for (cov in colnames(group.covs)) {
      scores.df[st:en, cov] <- group.covs[i,cov]
    }
    
    st <- st + nrow(score.matrix)
    en <- en + nrow(score.matrix)
  }
  
  if (add.direction) {
    scores.df$Direction <- sapply(strsplit(scores.df$Group, ';'), '[', 2)
    scores.df$Group <- sapply(strsplit(scores.df$Group, ';'), '[', 1)
  }
  
  return(scores.df)
}

flatten.score.matrix <- cmpfun(flatten.score.matrix)


## Helper function for p-value calculation.
get.multi.bins <- function(scores.df, all.covs, n.bins, use.quantiles = F, bin.direction = F, bin.group = F) {
  if (length(n.bins) == 1) {
    n.bins <- rep(n.bins, length(all.covs))
  }
  names(n.bins) <- all.covs
  
  bin.cov.list <- c()
  for (cov in all.covs) {
    bin.cov <- paste(cov, 'binned', sep = '_')
    if (use.quantiles) {
      bin.cov.list[[bin.cov]] <- factor(quantileCut(scores.df[[cov]], n = n.bins[[cov]], labels = F, include.lowest = T))
    } else {
      bin.cov.list[[bin.cov]] <- factor(cut(scores.df[[cov]], breaks = n.bins[[cov]], labels = F, include.lowest = T))
    }
  }
  
  if (bin.direction) {
    bin.cov.list[['Direction']] <- factor(scores.df$Direction)
  }
  
  if (bin.group) {
    bin.cov.list[['Group']] <- factor(scores.df$Group)
  }
  
  bin <- interaction(bin.cov.list, drop = T, sep = '_')
  bin.names <- levels(bin)
  
  bin.idx <- 1:length(bin.names)
  names(bin.idx) <- bin.names
  
  bin <- as.character(bin)
  
  scores.df$bin.index <- vapply(bin, function(x) bin.idx[[x]], 1)
  return(scores.df)
}

get.multi.bins <- cmpfun(get.multi.bins)


## Helper function for p-value calculation.
get.binned.list <- function(scores.df, output.name) {
  bin.names <- sort(unique(scores.df$bin.index))
  
  binned.dat <- lapply(bin.names, function(x) {
    scores <- scores.df[scores.df$bin.index == x,];
    if (nrow(scores) > 0) {
      return(sort(omitNaCpp(scores[[output.name]]), decreasing = F))
    } else {
      return(c())
    }
  })
  return(binned.dat)
}

get.binned.list <- cmpfun(get.binned.list)


## Helper function for p-value calculation.
fast.pvals <- function(r, binned.dat.up, binned.dat.dn) {
  bin.idx <- r[[1]]
  x <- r[[2]]
  if (x >= 0) {
    v <- binned.dat.up[[bin.idx]]
    p.val <- calcPvalGreaterCpp(v,x)
  } else {
    v <- binned.dat.dn[[bin.idx]]
    p.val <- calcPvalLessCpp(v,x)
  }
  return(p.val)
}

fast.pvals <- cmpfun(fast.pvals)


## Helper function for p-value calculation.
calc.emp.pvals <- function(scores.df, binned.dat, output.name, n.cores = 1, direction = c("both", "lower"), 
                           correct = "BH") {
  scores.mat <- as.matrix(scores.df[c('bin.index', output.name)])
  
  binned.dat.up <- lapply(binned.dat, function(v) v[v >= 0])
  binned.dat.dn <- lapply(binned.dat, function(v) v[v <= 0])
  rm(binned.dat)
  
  if (n.cores > 1) {
    cl <- makeCluster(n.cores, type = 'SOCK')
    clusterExport(cl, c('binned.dat.up', 'binned.dat.dn'), envir = environment())
    source.log <- parLapply(cl, 1:n.cores, function(x) {
      require('Rcpp')
      sourceCpp('hESC_cpp_functions.cpp');
    })
    print('Starting P-Value Calculations')
    coefs.pvals <- parApply(cl, scores.mat, 1, fast.pvals, binned.dat.up = binned.dat.up,
                            binned.dat.dn = binned.dat.dn)
    stopCluster(cl)
  } else {
    print('Starting P-Value Calculations')
    coefs.pvals <- apply(scores.mat, 1, fast.pvals, binned.dat.up = binned.dat.up,
                         binned.dat.dn = binned.dat.dn)
  }
  
  scores.df$p_val <- coefs.pvals
  if (correct == "qvalue") {
    scores.df$FDR <- qvalue(coefs.pvals)$qvalues
  } else {
    scores.df$FDR <- p.adjust(coefs.pvals, method = "BH")
  }
  scores.df <- scores.df[order(scores.df$p_val),]
  return(scores.df)
}

calc.emp.pvals <- cmpfun(calc.emp.pvals)