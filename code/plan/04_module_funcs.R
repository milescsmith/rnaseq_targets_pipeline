create_module_list <- function(module_file){
  module_tbl <- read_csv(module_file)

  module_list <-
    nest(
      .data = module_tbl,
      data = value
    ) %>%
    deframe() %>%
    map(pull, value)

  module_list
}

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
  ...,
  module_annotation
){
  module_list <- list(...)

  module_tbl <-
    map_dfr(
      module_list,
      function(x){
        enframe(x) %>% unnest(cols = "value")
        }
    ) %>%
    rename(
        module = name,
        gene   = value
    ) %>%
    inner_join(module_annotation)

  module_tbl
}

extract_module_genes <- function(
  module_table,
  exprs_mat,
  module_annotation = "Interferon"
  ){
    dplyr::group_by(
      .data = module_table,
      gene
    ) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::filter(
      type == module_annotation,
      gene %in% rownames(exprs_mat)
      ) %>%
    dplyr::select(module, gene) %>%
    arrange(module) %>%
    column_to_rownames("gene")
  }
