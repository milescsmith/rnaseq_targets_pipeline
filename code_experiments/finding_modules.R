library(reticulate)
use_condaenv('rtools', required = TRUE)
library(tidyverse)
#source("R/packages.R")
drake::loadd(vsd_exprs)

skd <- import("sklearn.decomposition")

top_vars <- matrixStats::rowSds(vsd_exprs) %>%
  `names<-`(rownames(vsd_exprs)) %>%
  enframe() %>%
  top_n(7500, value) %>%
  pull(name)

vsd_top <- vsd_exprs[top_vars,] %>% t()

standardize <- function(mat){
  (mat - mean(mat))/sd(mat)
}

ica = skd$FastICA(n_components=200L, random_state=NULL)

source = ica$fit_transform(t(standardize(vsd_top)))

ica_modules <- source %>%
  `rownames<-`(top_vars) %>%
  `colnames<-`(paste0("IC",1:200))

extract_modules <- function(mat, fdr_cutoff){
  map(1:ncol(mat), function(i){
    fdr_vals <- fdrtool::fdrtool(x = mat[,i],
                     plot = FALSE,
                     cutoff.method = "fndr")
    names(fdr_vals$qval[fdr_vals$qval < fdr_cutoff])
  })
  
}

ica_modules <- extract_modules(ica_modules, 1e-3) %>%
  `names<-`(paste0("IC", 1:length(ica_modules)))

ica_scores <- moduleScoreR::score_matrix(object = vsd_exprs,
                                         module_list = ica_modules) %>%
  as.matrix() %>%
  `rownames<-`(colnames(vsd_exprs)) %>%
  `colnames<-`(names(ica_modules))



loadd(annotation_info, group_pal)

library(WGCNA)
allowWGCNAThreads()

powers = 


plotDendroAndColors(dendro = wgcna_modules$dendrograms[[1]],
                    colors = wgcna_modules$colors[wgcna_modules$blockGenes[[1]]],
                    groupLabels = "Module colors",
                    main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

plotDendroAndColors(dendro = wgcna_modules$dendrograms[[2]],
                    colors = wgcna_modules$colors[wgcna_modules$blockGenes[[2]]], 
                    groupLabels = "Module colors",
                    main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

moduleColors = tibble(genes = colnames(vsd_top),
                      colors = wgcna_modules$colors)

plot(hclust(as.dist(1-cor(wgcna_modules$MEs)), method = "average"))
abline(h=0.25, col="red")

mergedModules = mergeCloseModules(exprData = vsd_top, 
                          colors = wgcna_modules$colors,
                          cutHeight = 0.3,
                          verbose = 3)

mergedColors = tibble(genes = colnames(vsd_top),
                      colors = mergedModules$colors)

mergedMEs = mergedModules$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(dendro = wgcna_modules$dendrograms[[1]],
                    colors = cbind(wgcna_modules$colors[wgcna_modules$blockGenes[[1]]],
                                   mergedModules$colors),
                    groupLabels = c("Initial modules","Merged modules"),
                    main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)


wgcna_module_genes <- wgcna_modules$colors %>%
  enframe(name = "gene",
          value = "module")

range_bound = 4
zscore_range <- c(-range_bound:range_bound)
pheatmap(mat = MEs[annotation_info %>%
                     as_tibble(rownames="sample_name") %>%
                     arrange(cluster) %>%
                     pull(sample_name),],
         scale = "column",cluster_rows = FALSE,
         annotation_row = annotation_info,
         annotation_colors = group_pal,
         gaps_row = annotation_info %>%
           group_by(cluster) %>%
           tally() %>%
           pull(n) %>%
           accumulate(sum) %>%
           `[`(1:length(.)-1),
         breaks = zscore_range,
         color = viridis(n = length(zscore_range)-1, option = "D"))

pheatmap(mat = MEs,
         scale = "column",
         cluster_rows = TRUE,
         annotation_row = annotation_info,
         annotation_colors = group_pal,
         breaks = zscore_range,
         color = viridis(n = length(zscore_range)-1, option = "D"))

pheatmap(mat = vsd_exprs[wgcna_module_genes %>%
                           filter(module == "turquoise") %>%
                           pull(gene),] %>%
           t(),
         scale = "column",
         annotation_row = annotation_info,
         annotation_colors = group_pal,
         gaps_row = annotation_info %>%
           group_by(cluster) %>%
           tally() %>%
           pull(n) %>%
           accumulate(sum) %>%
           `[`(1:length(.)-1),
         breaks = zscore_range,
         color = viridis(n = length(zscore_range)-1, option = "D"))


pheatmap(mat = ica_scores[annotation_info %>%
                     as_tibble(rownames="sample_name") %>%
                     arrange(cluster) %>%
                     pull(sample_name),],
         scale = "column",
         cluster_rows = FALSE,
         annotation_row = annotation_info,
         annotation_colors = group_pal,
         breaks = zscore_range,cutree_rows = 7,
         gaps_row = annotation_info %>%
           group_by(cluster) %>%
           tally() %>%
           pull(n) %>%
           accumulate(sum) %>%
           `[`(1:length(.)-1),
         color = viridis(n = length(zscore_range)-1, option = "D"))
