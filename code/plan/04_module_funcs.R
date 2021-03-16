extract_module_scores <- function(
  dds,
  ...
){
  module_names <- dots_list(...)

  colData(dds) %>%
    as_tibble(rownames = "sample_name") %>%
    select(
      sample_name,
      one_of(!!!module_names)
      )
}

create_module_table <- function(
  ...
){
  module_list <- dots_list(...)
  
  enframe(!!!module_list) %>%
    unnest(cols = "value") %>%
    rename(
      module = name,
      gene   = value
      ) %>%
    inner_join(module_annotation)
}

extract_module_genes <- function(
  module_table,
  module_annotation = "Interferon"
  ){
    dplyr::group_by(
      module_tbl,
      gene
    ) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::filter(
      type == "Interferon",
      gene %in% rownames(vsd_exprs)
      ) %>%
    dplyr::select(module, gene) %>%
    arrange(module) %>%
    column_to_rownames("gene")
  }