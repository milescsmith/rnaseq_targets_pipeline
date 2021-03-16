#### WGCNA ####
top_variable_genes <- function(exprs, n = 20000){
  top_vars =
    rowSds(exprs) %>%
    set_names(rownames(exprs)) %>%
    enframe() %>%
    top_n(n, value) %>%
    pull(name)

  vsd_top = exprs[top_vars,] %>%
    t()

  vsd_top_float <- `storage.mode<-`(vsd_top, "numeric")
  vsd_top_float
}


module_gsea <- function(module_genes, module_of_interest){
  enriched_module_genes <-
    module_genes %>%
      filter(module == module_of_interest) %>%
      pull(hugo) %>%
      enricher(gene = .,
                TERM2GENE = c5)
  
  if(!is.null(enriched_module_genes)){
    mutate(
      .data = slot(enriched_module_genes,'result'),
      module = module_of_interest
      )
  }
  
  enriched_module_genes
}


module_gsea_plots <- function(enriched_genes){
  filter(
    .data = enriched_genes,
    p.adjust < 0.05,
    module != "grey") %>%
  mutate(
    GeneRatio = map_dbl(
      .x = GeneRatio,
      .f = function(i){
        j = str_split(i, "/") %>%
          magrittr::extract2(1) %>%
          as.double()
        j[[1]]/j[[2]]
      }),
    ID =
      str_replace_all(
        string = ID,
        pattern = "_",
        replacement = " "
        ),
    module = paste0("ME", module)
  ) %>%
  group_by(module) %>%
  top_n(
    n = 5,
    wt = GeneRatio
    ) %>%
  sample_n(
    size = 5,
    replace = TRUE
    ) %>%
  distinct() %>%
  ungroup() %>%
  arrange(
    module,
    GeneRatio
    ) %>%
  mutate(order = row_number())
}
