library(igraph)
library(perturbLM)

## Format coefficient matrix from all media conditions
coefs.df <- read.table("up-tf-all_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.endo <- read.table("up-tf-endo_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)
coefs.df.multi <- read.table("up-tf-multi_regression_trimmed.pvals.tsv", sep = "\t", header = T, stringsAsFactors = F)

genotype.summary <- read.table("up-tf-all_genotype_summary.tsv", sep = "\t")
genotype.summary.endo <- read.table("up-tf-endo_genotype_summary.tsv", sep = "\t")
genotype.summary.multi <- read.table("up-tf-multi_genotype_summary.tsv", sep = "\t")

sig.genotypes <- union(rownames(subset(genotype.summary, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                            rownames(subset(genotype.summary, n_diff_genes >= 100)))
sig.genotypes.endo <- union(rownames(subset(genotype.summary.endo, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                            rownames(subset(genotype.summary.endo, n_diff_genes >= 100)))
sig.genotypes.multi <- union(rownames(subset(genotype.summary.multi, chisq.fdr < 1e-6 & min.cluster.fdr < 1e-6)),
                            rownames(subset(genotype.summary.multi, n_diff_genes >= 100)))

coefs.matrix.stem <- UnflattenDataframe(coefs.df, "cf", row.col = "Gene", col.col = "Group")
coefs.matrix.stem <- coefs.matrix.stem[, sig.genotypes]

coefs.matrix.endo <- UnflattenDataframe(coefs.df.endo, "cf", row.col = "Gene", col.col = "Group")
coefs.matrix.endo <- coefs.matrix.endo[, sig.genotypes.endo]

coefs.matrix.multi <- UnflattenDataframe(coefs.df.multi, "cf", row.col = "Gene", col.col = "Group")
coefs.matrix.multi <- coefs.matrix.multi[, sig.genotypes.multi]

## Create gene-gene network
gene.coefs.snn <- BuildSNN(coefs.matrix.stem, k.param = 30, k.scale = 10, prune.SNN = 1/30)
gene_func_graph <- graph_from_adjacency_matrix(as.matrix(gene.coefs.snn), mode = "undirected", weighted = T)

## Segment genes into modules
gene_func_clusters <- cluster_louvain(gene_func_graph, weights = E(gene_func_graph)$weights)
gene_clusters <- sapply(membership(gene_func_clusters), function(i) paste("GM", i, sep = ""))

gene_clusters_list <- lapply(unique(gene_clusters), function(i) names(gene_clusters[gene_clusters == i]))
names(gene_clusters_list) <- unique(gene_clusters)
print(sapply(gene_clusters_list[order(names(gene_clusters_list))], length))

# write.genesets(gene_clusters_list, "up-tf_functional_gene_modules.gmt")
# write_graph(gene_func_graph, "up-tf_functional_gene_network.edgelist", "ncol")

## Load the SNN network and gene modules generated for the paper
gene_func_graph <- read.graph("up-tf_functional_gene_network.edgelist", format="ncol")
gene_clusters_list <- LoadGenesets("up-tf_functional_gene_modules.gmt")
gene_clusters <- FlattenGenotypeList(gene_clusters_list)

## GO enrichment
genesets <- LoadGenesets("go_bp_genesets.gmt")

bg.genes <- unique(unlist(gene_clusters_list, F, F))
genesets <- FilterGenesets(genesets, bg.genes, min.size = 30, max.size = 500)
names(genesets) <- sapply(names(genesets), function(x) gsub("GO_", "", x))

gene.clusters.enrich <- MultipleFisherEnrich(gene_clusters_list, genesets)
top.gene.clusters.enrich <- subset(gene.clusters.enrich, p.val < 0.05)
top.gene.clusters.enrich <- Reduce(rbind, by(top.gene.clusters.enrich, top.gene.clusters.enrich$Group, head, n = 3))

max.score <- 30
heat.df <- subset(gene.clusters.enrich, genesets %in% top.gene.clusters.enrich$genesets)
heat.df$Score <- -log(heat.df$p.val)
heat.mat <- UnflattenDataframe(heat.df, output.name = "Score", row.col = 'genesets', col.col = 'Group')
heat.mat[abs(heat.mat) > max.score] <- max.score
heat.mat <- heat.mat[unique(top.gene.clusters.enrich$genesets),]

# pdf("gene_module_go_enrichment.pdf", width = 11, height = 6)
ggheat(heat.mat, clustering = "none", labRow = T, x.lab.size = 12, y.lab.size = 10)
# dev.off()

# write.table(gene.clusters.enrich, file = "up-tf_gene_module_enrichment.tsv", sep = "\t")

## Identify which genes in each module are part of enriched genesets
cluster.marker.genes <- lapply(names(gene_clusters_list), function(module) {
  module.genes <- gene_clusters_list[[module]]
  module.genesets <- subset(top.gene.clusters.enrich, Group == module)$genesets
  
  module.genes <- module.genes[module.genes %in% unlist(genesets[module.genesets], F, F)]
  genes2sets <-lapply(module.genes, function(gene) { 
    sets <- sapply(genesets[module.genesets], function(set) gene %in% set)
    set.names <- names(sets[sets == T])
    set.names
  })
  genes2sets.counts <- sapply(genes2sets, length)
  genes2sets <- sapply(genes2sets, function(x) paste(x, collapse = ","))
  
  res.df <- data.frame(Module = module, Gene = module.genes, Genesets = genes2sets, N_Sets = genes2sets.counts)
  res.df[order(res.df$N_Sets, decreasing = T),]
})
cluster.marker.genes <- do.call("rbind", cluster.marker.genes)
rownames(cluster.marker.genes) <- NULL
cluster.marker.genes <- cluster.marker.genes[order(cluster.marker.genes$Module),]
gene_membership <- data.frame(gene = names(gene_clusters), module = gene_clusters)

# write.table(gene_membership, "gene_module_membership.tsv", sep = "\t", row.names = F)
# write.table(cluster.marker.genes, file = "up-tf_gene_module_geneset_assignments.tsv", sep = "\t")

## Visualize gene module relationships
gene_clusters_integer <- sapply(gene_clusters, function(x) gsub("GM", "", x))
V(gene_func_graph)$size <- 1
gene_clusters_graph <- contract.vertices(gene_func_graph, gene_clusters_integer[names(V(gene_func_graph))],
                                         vertex.attr.comb = list(size = "sum", "ignore"))
gene_clusters_graph <- simplify(gene_clusters_graph, remove.loops = T, edge.attr.comb = list(weight = "sum", "ignore"))
E(gene_clusters_graph)$weight <- E(gene_clusters_graph)$weight/mean(E(gene_clusters_graph)$weight)
V(gene_clusters_graph)$size <- V(gene_clusters_graph)$size/mean(V(gene_clusters_graph)$size)
gene_clusters_graph <- delete_edges(gene_clusters_graph, which(E(gene_clusters_graph)$weight < 0.1))
l <- layout_with_fr(gene_clusters_graph, weight = E(gene_clusters_graph)$weight*0.25)

# pdf("gene_module_network_labels.pdf", width = 5.5, height = 5.5)
plot(gene_clusters_graph, layout = l, vertex.size = (V(gene_clusters_graph)$size)^0.5*15, edge.width = E(gene_clusters_graph)$weight,
     vertex.label.cex = 0.75, vertex.color = "slategrey")
# dev.off()

## Analyze TF effects on gene modules
coefs.matrices.list <- list(Stem = coefs.matrix.stem, Endo = coefs.matrix.endo, Multi = coefs.matrix.multi)
gene_module_mapping <- read.table("up-tf_gene_module_mapping.txt", sep = "\t", header = T)

tfs.modules.matrices.list <- lapply(coefs.matrices.list, function(X) {
  CalcGeneModuleEffect(X, gene_clusters_list, gene_module_mapping, min.coef = 0.025)
})
tfs.modules.matrix <- do.call(cbind, tfs.modules.matrices.list)

## Normalize to mCherry
max.score <- 0.075
min.score <- -0.075
tfs.modules.matrix[tfs.modules.matrix > max.score] <- max.score
tfs.modules.matrix[tfs.modules.matrix < min.score] <- min.score

# pdf("tf_gene_module_effect_norm.pdf", width = 15, height = 4.5)
ggheat(tfs.modules.matrix, clustering = "row", labRow = T, x.lab.size = 13.5, y.lab.size = 13.5)
# dev.off()

## Effect of specific TFs on pluripotency subnetwork
TF <- "KLF4"
#TF <- "SNAI2"
seed.nodes <- c("SOX2", "POU5F1", "NANOG") ## Pluripotency subnetwork
core.nodes <- c(seed.nodes, "DNMT3B", "DPPA4", "SALL2")

core_subgraph <- induced_subgraph(gene_func_graph, vids = core.nodes)
core_subgraph <- simplify(core_subgraph)
core_subgraph <- delete_vertices(core_subgraph, degree(core_subgraph) == 0)
l <- layout_with_fr(core_subgraph, weights = E(core_subgraph)$weight*10)

tfs.coefs <- coefs.matrix.stem[,TF]
V(core_subgraph)$size <- abs(tfs.coefs[names(V(core_subgraph))])
V(core_subgraph)$color <- sapply(tfs.coefs[names(V(core_subgraph))], function(x) ifelse(sign(x) > 0, "tomato", "skyblue"))

# pdf("pluripotency_subnetwork.pdf", width = 5, height = 5)
plot(core_subgraph, vertex.size = V(core_subgraph)$size*100, vertex.color = V(core_subgraph)$color,
     layout = l, edge.width = E(core_subgraph)$weight * 10)
# dev.off()