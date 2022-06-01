library(biomaRt)

mart <-
  useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl"
  )

up_translation_table <-
  getBM(
    attributes = c("hgnc_symbol", "entrezgene_id"),
    filters    = "hgnc_symbol",
    values     = upregulated$gene,
    mart       = mart
    )


ggo <- groupGO(gene = as.character(upreg_entrez_id), OrgDb = "org.Hs.eg.db", ont = "BP", level = 3, readable = TRUE)

upreg_entrez_id <- up_translation_table %>% filter(!is.na(entrezgene_id)) %>% pull(entrezgene_id)
up_ego_bp <- enrichGO(
  gene = as.character(upreg_entrez_id),
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  readable = TRUE
  )

goplot(up_ego_bp)

up_ego_mf <- enrichGO(
  gene = as.character(upreg_entrez_id),
  OrgDb = "org.Hs.eg.db",
  ont = "MF",
  readable = TRUE
)

goplot(up_ego_mf)

down_translation_table <-
  getBM(
    attributes = c("hgnc_symbol", "entrezgene_id"),
    filters    = "hgnc_symbol",
    values     = downregulated$gene,
    mart       = mart
  )

downreg_entrez_id <- down_translation_table %>% filter(!is.na(entrezgene_id)) %>% pull(entrezgene_id)
down_ego_bp <- enrichGO(
  gene = as.character(downreg_entrez_id),
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  readable = TRUE
)

goplot(down_ego_bp)

down_ego_mf <- enrichGO(
  gene = as.character(downreg_entrez_id),
  OrgDb = "org.Hs.eg.db",
  ont = "MF",
  readable = TRUE
)

goplot(down_ego_mf)

hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)
genes <- upregulated$gene
mgeneSim(genes, semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)


up_era <- enrichPathway(
  gene = as.character(upreg_entrez_id),
  pvalueCutoff = 0.05,
  readable = TRUE
)

heatplot(up_era, foldChange = upregulated %>% dplyr::select(gene, log2FoldChange) %>% tibble::deframe())
cnetplot(up_era)

up_era <- enrichPathway(
  gene = as.character(upreg_entrez_id),
  pvalueCutoff = 0.05,
  readable = TRUE
)

upFoldChange <-
  upregulated %>%
  dplyr::select(gene, log2FoldChange) %>%
  tibble::deframe()

{
  png(file = "up_ra.png", width = 600, height = 600, units = "px")
  cnetplot(
    x = up_era,
    foldChange = upFoldChange,
    colorEdge = TRUE,
    layout = "fr",
    showCategory = 10,
    cex_label_gene = 0.5,
    color_gene = "#000000",
    shadowtext = "category"
  ) +
  guides(
    size = "none",
    edge_colour = "none"
  ) +
  labs(title = "Reactome pathways of genes upregulated in SARS2-CoV-2 infected")
  dev.off()
}

down_era <- enrichPathway(
  gene = as.character(downreg_entrez_id),
  pvalueCutoff = 0.05,
  readable = TRUE
)

downFoldChange <-
  downregulated %>%
  dplyr::select(gene, log2FoldChange) %>%
  tibble::deframe()

{
  png(file = "down_ra.png", width = 600, height = 600, units = "px")
  cnetplot(
    x = down_era,
    foldChange = downFoldChange,
    colorEdge = TRUE,
    layout = "fr",
    showCategory = 13,
    cex_label_gene = 0.75,
    color_gene = "#000000",
    shadowtext = "category"
  ) +
  guides(
    size = "none",
    edge_colour = "none"
  ) +
  labs(title = "Reactome pathways of genes downregulated in SARS2-CoV-2 infected")
  dev.off()
}

{
  png(file = "down_ra_emapplot.png", width = 600, height = 600, units = "px")
  emapplot(
    x = pairwise_termsim(x = down_era),
    showCategory = 10
    )
  dev.off()
}

{
  png(file = "up_ra_emapplot.png", width = 600, height = 600, units = "px")
  emapplot(
    x = pairwise_termsim(x = up_era),
    showCategory = 10
  )
  dev.off()
}

cnetplot(up_enrichment$`sle - control`$gene_ontology, foldChange = )
