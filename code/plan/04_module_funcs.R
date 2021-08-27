#' @title extract_module_scores
#'
#' @description Extract the module scores from the
#' the metadata associated with a data object
#'
#' @param object A data object (either DESeqDataSet or DGEList)
#' @param ... a list of module names
#'
#' @return A tibble with the sample names and module values
#' @export
#'
extract_module_scores <- function(object, ...){
    UseMethod("extract_module_scores", object)
  }

#' @rdname extract_module_scores
#' @method extract_module_scores DESeqDataSet
#' @importFrom purrr chuck
#' @importFrom rlang dots_list
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tidyselect one_of
#' @return tibble
#'
extract_module_scores.DESeqDataSet <-
  function(
    object,
    ...
    ){
    module_names <- rlang::dots_list(...)
    object %>%
      purrr::chuck("colData") %>%
      tibble::as_tibble(
        rownames = "sample_name"
        ) %>%
      dplyr::select(
        sample_name,
        tidyselect::one_of(!!!module_names)
      )
    }

#' @rdname extract_module_scores
#' @method extract_module_scores DGEList
#' @importFrom purrr chuck
#' @importFrom rlang dots_list
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom tidyselect one_of
#' @return tibble
#'
extract_module_scores.DGEList <-
  function(
    object,
    ...
  ){
    module_names <- rlang::dots_list(...)
    object %>%
      purrr::chuck("samples") %>%
      tibble::as_tibble(
        rownames = "sample_name"
      )
      dplyr::select(
        sample_name,
        tidyselect::one_of(!!!module_names)
      )
}

create_module_list <- function(module_file){
  readr::read_csv(module_file) %>%
    tidyr::nest(data = value) %>%
    tibble::deframe() %>%
    purrr::map(.f = dplyr::pull, value)
}


create_module_table <- function(
  ...,
  module_annotation
){
  module_list <- rlang::list(...)

  module_tbl <-
    purrr::map_dfr(
      module_list,
      \(x) tibble::enframe(x) %>% tidyr::unnest(cols = "value")
    ) %>%
    dplyr::rename(
      module = name,
      gene   = value
    ) %>%
    dplyr::inner_join(module_annotation)

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
    dplyr::ungroup() %>%
    dplyr::filter(
      type == module_annotation,
      gene %in% rownames(exprs_mat)
      ) %>%
    dplyr::select(module, gene) %>%
    dplyr::arrange(module) %>%
    tibble::column_to_rownames("gene")
}
