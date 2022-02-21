#### WGCNA ####
top_variable_genes <-
  function(
    exprs,
    n = 20000
    ){
    top_vars =
      matrixStats::rowSds(exprs) %>%
      rlang::set_names(rownames(exprs)) %>%
      tibble::enframe() %>%
      dplyr::top_n(n, value) %>%
      dplyr::pull(name)

      vsd_top =
        exprs[top_vars,] %>%
        t()

        # vsd_top_float <- `storage.mode<-`(vsd_top, "numeric")
        # vsd_top_float

      vsd_top
  }


module_gsea <-
  function(
    module_genes,
    module_of_interest,
    target_species = "org.Hs.eg.db"
  ){
    enriched_module_genes <-
      module_genes |>
      dplyr::filter(module == module_of_interest) |>
      dplyr::pull(hugo) |>
      clusterProfiler::enrichGO(
        OrgDb = target_species,
        keyType = "SYMBOL",
        ont = "ALL"
      )
    if(!is.null(enriched_module_genes)){
      enriched_module_genes <-
        dplyr::mutate(
          .data = enriched_module_genes@result,
          module = module_of_interest
        )
    }

    enriched_module_genes
  }


module_gsea_plots <-
  function(
    enriched_genes
  ){
    dplyr::filter(
      .data = enriched_genes,
      p.adjust < 0.05,
      module != "grey") %>%
      dplyr::mutate(
        GeneRatio = purrr::map_dbl(
          .x = GeneRatio,
          .f = function(i){
            j = stringr::str_split(i, "/") %>%
              magrittr::extract2(1) %>%
              as.double()
            j[[1]]/j[[2]]
          }),
        ID =
          stringr::str_replace_all(
            string = ID,
            pattern = "_",
            replacement = " "
          ),
        module = paste0("ME", module)
      ) %>%
      dplyr::group_by(module) %>%
      dplyr::top_n(
        n = 5,
        wt = GeneRatio
      ) %>%
      dplyr::sample_n(
        size = 5,
        replace = TRUE
      ) %>%
      dplyr::distinct() %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        module,
        GeneRatio
      ) %>%
      dplyr::mutate(order = dplyr::row_number())
  }

find_softPower <- function(sft){
  if (is.na(sft$powerEstimate)){
    scale_free_topo_fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
    powerEstimate <- which(scale_free_topo_fit == max(scale_free_topo_fit))
  } else {
    powerEstimate <- sft$powerEstimate
  }

  powerEstimate
}
