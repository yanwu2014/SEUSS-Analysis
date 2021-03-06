library(igraph)
library(perturbLM)
library(swne)

#### Load data and prep combined coefficient matrix ####

## Read in linear model coefficients
coefs.df.stem <- read.table("up-tf-stem_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.endo <- read.table("up-tf-endo_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.multi <- read.table("up-tf-multi_regression.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)

## Read in all summarized results
genotype.summary.stem <- read.table("up-tf-stem_summary.tsv", sep = "\t")
genotype.summary.endo <- read.table("up-tf-endo_summary.tsv", sep = "\t")
genotype.summary.multi <- read.table("up-tf-multi_summary.tsv", sep = "\t")

## Filter for significant genotypes
min.cl.pval <- 1e-12
min.diff.genes <- 40

ctrl.cl.stem <- genotype.summary.stem["mCherry", "min_cluster"]
sig.genotypes.stem <- rownames(subset(genotype.summary.stem, (abs(min_cluster_pval) < min.cl.pval & min_cluster != ctrl.cl.stem) 
                                      | n_diff_genes >= min.diff.genes))
ctrl.cl.endo <- genotype.summary.endo["mCherry", "min_cluster"]
sig.genotypes.endo <- rownames(subset(genotype.summary.endo, (abs(min_cluster_pval) < min.cl.pval & min_cluster != ctrl.cl.endo) 
                                      | n_diff_genes >= min.diff.genes))
ctrl.cl.multi <- genotype.summary.multi["mCherry", "min_cluster"]
sig.genotypes.multi <- rownames(subset(genotype.summary.multi, (abs(min_cluster_pval) < min.cl.pval & min_cluster != ctrl.cl.multi) 
                                       | n_diff_genes >= min.diff.genes))

## Re-create linear coefficient matrices
coefs.matrix.stem <- UnflattenDataframe(coefs.df.stem, "cf")
coefs.matrix.endo <- UnflattenDataframe(coefs.df.endo, "cf")
coefs.matrix.multi <- UnflattenDataframe(coefs.df.multi, "cf")

## Combine matrices and filter for genes that are significant in at least one genotype
min.abs.cf <- 0.025
max.fdr <- 0.2

genes.all <- Reduce(intersect, list(rownames(coefs.matrix.stem), rownames(coefs.matrix.endo), rownames(coefs.matrix.multi)))
genes.ix <- unique(c(subset(coefs.df.stem, abs(cf) > min.abs.cf & FDR < max.fdr & Group %in% sig.genotypes.stem)$Gene,
                      subset(coefs.df.endo, abs(cf) > min.abs.cf & FDR < max.fdr & Group %in% sig.genotypes.endo)$Gene))
genes.ix <- intersect(genes.ix, genes.all)

coefs.matrix <- cbind(coefs.matrix.stem[genes.ix, sig.genotypes.stem],
                      coefs.matrix.endo[genes.ix, sig.genotypes.endo])
dim(coefs.matrix)

#### Create gene SNN graph and segment genes into modules ####

## Create gene-gene network
gene.coefs.snn <- CalcSNN(t(coefs.matrix), k = 25, prune.SNN = 1/10)
rownames(gene.coefs.snn) <- colnames(gene.coefs.snn) <- rownames(coefs.matrix)
gene.coefs.graph <- graph_from_adjacency_matrix(as.matrix(gene.coefs.snn), mode = "undirected", weighted = T)
V(gene.coefs.graph)$name <- rownames(gene.coefs.snn)

## Segment genes into modules
gene.modules <- cluster_louvain(gene.coefs.graph, weights = E(gene.coefs.graph)$weights)
gene.modules <- membership(gene.modules)

gene.modules.counts <- table(gene.modules)
gene.modules.counts <- gene.modules.counts[gene.modules.counts > 20]
gene.modules <- factor(gene.modules[gene.modules %in% names(gene.modules.counts)])
levels(gene.modules) <- paste("GM", 1:nlevels(gene.modules), sep = "")
table(gene.modules)

gene.coefs.graph <- induced_subgraph(gene.coefs.graph, vids = which(V(gene.coefs.graph)$name %in% names(gene.modules)))
V(gene.coefs.graph)$size <- 1
write_graph(gene.coefs.graph, "up-tf_gene_network.edgelist", "ncol")

gene.modules.list <- UnflattenGroups(gene.modules)
WriteGenesets(gene.modules.list, "up-tf_gene_modules.gmt")

#### Analyze genotype effects on gene modules ####

## Load the SNN network and gene modules generated for the paper
gene.coefs.graph <- read_graph("up-tf_gene_network.edgelist", format = "ncol")
gene.modules.list <- LoadGenesets("up-tf_gene_modules.gmt")
gene.modules <- FlattenGenotypeList(gene.modules.list)
table(gene.modules)

module.names <- read.table("up-tf_module_names.txt", sep = "\t", header = T)
colnames(module.names) <- c("Module", "Description")
# module.names <- NULL

coefs.matrix.list <- list(coefs.matrix.stem[genes.ix, unique(c(sig.genotypes.stem, "mCherry-int-ctrl"))],
                          coefs.matrix.endo[genes.ix, unique(c(sig.genotypes.endo, "mCherry-int-ctrl"))],
                          coefs.matrix.multi[genes.ix, unique(c(sig.genotypes.multi, "mCherry-int-ctrl"))])
tfs.modules.matrix <- lapply(coefs.matrix.list, function(X) {
  CalcGeneModuleEffect(X, gene.modules.list, module.names, min.coef = 0.01)
})
tfs.modules.matrix <- do.call(cbind, tfs.modules.matrix)
tfs.modules.matrix <- tfs.modules.matrix[apply(abs(tfs.modules.matrix), 1, max) > 0.05,]

max.score <- 0.1
tfs.modules.matrix[tfs.modules.matrix > max.score] <- max.score
tfs.modules.matrix[tfs.modules.matrix < -1*max.score] <- -1*max.score

# pdf("up-tf_gene_module_effect.pdf", width = 12, height = 4.0)
ggHeat(tfs.modules.matrix, clustering = "row", labRow = T, x.lab.size = 14, y.lab.size = 14)
# dev.off()


#### Interpret gene modules ####

## Load genesets
genesets <- LoadGenesets("hallmark_go_bp_genesets.gmt")
genesets <- FilterGenesets(genesets, genes.all, min.size = 15, max.size = 1000)
names(genesets) <- sapply(names(genesets), function(x) gsub("GO_", "", x))
length(genesets)

## Run enrichment
modules.enrich <- MultipleFisherEnrich(gene.modules.list, genesets)
top.modules.enrich <- subset(modules.enrich, p.val < 0.05)
top.modules.enrich <- Reduce(rbind, by(top.modules.enrich, top.modules.enrich$Group, head, n = 10))
top.modules.enrich <- top.modules.enrich[order(as.character(top.modules.enrich$Group)),]
table(top.modules.enrich$Group)

write.table(top.modules.enrich, file = "up-tf_module_geneset_enrich.tsv", sep = "\t")

max.lp <- 20
top.modules.enrich <- subset(top.modules.enrich, Group %in% rownames(tfs.modules.matrix))
top.modules.enrich <- Reduce(rbind, by(top.modules.enrich, top.modules.enrich$Group, head, n = 4))
heat.df <- subset(modules.enrich, genesets %in% top.modules.enrich$genesets)
heat.df$lp <- -log(heat.df$p.val)
heat.mat <- UnflattenDataframe(heat.df, output.name = "lp", row.col = 'genesets', col.col = 'Group')
heat.mat[abs(heat.mat) > max.lp] <- max.lp
heat.mat <- heat.mat[unique(top.modules.enrich$genesets), rownames(tfs.modules.matrix)]

pdf("up-tf_module_geneset_enrich_heatmap.pdf", width = 10, height = 8)
ggHeat(heat.mat, clustering = "both", labRow = T, x.lab.size = 12, y.lab.size = 10)
dev.off()

## Identify which genes in each module are part of enriched genesets
top.modules.enrich <- subset(modules.enrich, p.val < 0.05 & Group %in% rownames(tfs.modules.matrix))
module.genes <- lapply(rownames(tfs.modules.matrix), function(m) {
  m.enrich <- subset(top.modules.enrich, Group == m)
  m.enrich <- m.enrich[order(m.enrich$p.val),]
  m.genesets <- m.enrich$genesets
  m.genesets.p <- m.enrich$p.val; names(m.genesets.p) <- m.genesets;
  
  m.genes <- gene.modules.list[[m]]
  m.genes <- m.genes[m.genes %in% unlist(genesets[m.genesets], F, F)]
  genes.sets <- lapply(m.genes, function(gene) {
    sets <- sapply(genesets[m.genesets], function(set) gene %in% set)
    sets <- names(sets[sets == T])
    m.genesets.p[sets]
  })
  genes.sets.counts <- sapply(genes.sets, length)
  genes.sets.lp <- sapply(genes.sets, function(p) sum(-log(p)))
  genes.sets.names <- sapply(genes.sets, function(x) paste(names(x), collapse = ","))
  
  df <- data.frame(module = m, gene = m.genes, lp_sum = genes.sets.lp, n_genesets = genes.sets.counts, 
                   genesets = genes.sets.names)
  df[order(df$lp_sum, decreasing = T),]
})
module.genes <- do.call("rbind", module.genes); rownames(module.genes) <- NULL;
module.genes <- Reduce(rbind, by(module.genes, module.genes$module, head, n = 20))
module.genes <- module.genes[order(as.character(module.genes$module)),]

write.table(module.genes, file = "up-tf_module_geneset_genes.tsv", sep = "\t")


## Look at top genes in top genotypes for each gene module
module.coef.genes <- lapply(rownames(tfs.modules.matrix), function(m) {
  m.genes <- gene.modules.list[[m]]
  abs.coef <- rowMeans(abs(coefs.matrix[m.genes,]))
  df <- data.frame(module = m, gene = m.genes, abs_coef = abs.coef)
  df[order(df$abs_coef, decreasing = T),]
})
module.coef.genes <- do.call("rbind", module.coef.genes); rownames(module.coef.genes) <- NULL;
module.coef.genes <- Reduce(rbind, by(module.coef.genes, module.coef.genes$module, head, n = 20))
module.coef.genes <- module.coef.genes[order(as.character(module.coef.genes$module)),]

write.table(module.coef.genes, file = "up-tf_module_coef_genes.tsv", sep = "\t")



#### Plot gene module relationships ####
gene.modules <- gene.modules[gene.modules %in% module.names$Module]
gene.coefs.graph <- induced_subgraph(gene.coefs.graph, vids = which(V(gene.coefs.graph)$name %in% names(gene.modules)))
V(gene.coefs.graph)$size <- 1
gene.modules.int <- sapply(as.character(gene.modules), function(x) gsub("GM", "", x))
names(gene.modules.int) <- names(gene.modules)

modules.graph <- contract(gene.coefs.graph, gene.modules.int[V(gene.coefs.graph)$name],
                          vertex.attr.comb = list(size = "sum", "ignore"))
modules.graph <- simplify(modules.graph, remove.loops = T, edge.attr.comb = list(weight = "sum", "ignore"))
modules.graph <- delete.vertices(simplify(modules.graph), degree(modules.graph) == 0)
V(modules.graph)$label <- as.character(module.names$Description)

E(modules.graph)$weight <- E(modules.graph)$weight/mean(E(modules.graph)$weight)
V(modules.graph)$size <- V(modules.graph)$size/mean(V(modules.graph)$size)

modules.graph <- delete_edges(modules.graph, which(E(modules.graph)$weight < 0.05))
modules.graph.layout <- layout_with_fr(modules.graph, weight = E(modules.graph)$weight*0.5)

library(RColorBrewer)
pdf("up-tf_module_network_labels.pdf", width = 6, height = 6)
plot(modules.graph, layout = modules.graph.layout, vertex.size = sqrt(V(modules.graph)$size)*20, 
     edge.width = E(modules.graph)$weight, vertex.label.cex = 1, vertex.label = V(modules.graph)$label,
     vertex.color = brewer.pal(gorder(modules.graph), "Spectral"))
dev.off()

pdf("up-tf_module_network_nolabels.pdf", width = 6, height = 6)
plot(modules.graph, layout = modules.graph.layout, vertex.size = sqrt(V(modules.graph)$size)*25, 
     edge.width = E(modules.graph)$weight, vertex.label = NA, 
     vertex.color = brewer.pal(gorder(modules.graph), "Spectral"))
dev.off()

