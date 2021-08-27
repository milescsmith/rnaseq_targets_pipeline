#' @title create_results_list
#'
#' @param object
#' @param ...
#'
#' @return
#' @export
create_results_list <- function(object, ...){
  UseMethod("create_results_list")
}

#' @rdname create_results_list
#' @method create_results_list DESeqDataSet
#' @importFrom DESeq2 lfcShrink
#' @importFrom rlang set_names
#' @importFrom purrr map map_chr
#' @importFrom dplyr mutate inner_join filter pull
#' @return list
create_results_list.DESeqDataSet <- function(
  comparison_list,
  object,
  comparison_grouping_variable,
  design = NULL
  ){
  purrr::map(comparison_list, function(i) {
    DESeq2::lfcShrink(
      dds = object,
      coef = i,
      parallel = TRUE,
      type = "apeglm"
      ) %>%
    tibble::as_tibble(rownames="gene")
  }) %>%
    rlang::set_names(
      purrr::map_chr(
        comparison_list,
        str_remove,
        pattern = paste0(comparison_grouping_variable, "_")
      )
  )
}

# TODO: FINISH HIM!
#' @rdname create_results_list
#' @method create_results_list DGEList
#' @importFrom magrittr use_series
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest getCounts
#' @importFrom rlang set_names
#' @importFrom purrr map map_chr chuck
#' @importFrom tibble as_tibble
#' @importFrom matrixStats rowMeans2
#' @importFrom rstatix adjust_pvalue
#' @importFrom dplyr mutate inner_join filter pull rename relocate
#' @return list
create_results_list.DGEList <- function(
  comparison_list,
  object,
  comparison_grouping_variable,
  design = NULL
  ){

  if (is.null(magrittr::use_series(object, common.dispersion))){
    message("No dispersion values found.")
    if (is.null(design)){
      stop("No design formula or matrix was provided.  The design is necessary
           to calculate dispersion")
    } else {
      message("Calculating dispersion")
      object <-
        edgeR::estimateDisp(
          y = object,
          design = design,
          robust = TRUE
        )
    }
  }

  comparisons <-
    map(
      comparison_list,
      \(x) str_split(x, pattern =  " - ") |>
        chuck(1) |>
        chuck(1) %>%
        paste0(comparison_grouping_variable, .)
    )

  purrr::map(
    .x = comparisons,
    .f = \(i) {
    edgeR::glmQLFit(
      y           = object,
      design      = mm,
      prior.count = 3
      ) %>%
    edgeR::glmQLFTest(coef = i) %>%
    magrittr::use_series(table) %>%
    tibble::as_tibble(rownames="gene") %>%
    dplyr::rename(
      log2FoldChange = logFC,
      pvalue = PValue
    ) %>%
    rstatix::adjust_pvalue(
      p.col      = "pvalue",
      output.col = "padj",
      method     = "fdr"
    ) %>%
    dplyr::mutate(baseMean = matrixStats::rowMeans2(edgeR::getCounts(object))) %>%
    dplyr::relocate(
      gene,
      baseMean,
      log2FoldChange,
      F,
      pvalue,
      padj
      )
  }) %>%
    rlang::set_names(
      purrr::map_chr(
        comparison_list,
        str_remove,
        pattern = paste0(comparison_grouping_variable, "_")
      )
  )
}


create_deg_tables <- function(
  deg_res,
  comparison_list,
  grouping_variable,
  direction=c("up","down")
){
  direction <- match.arg(arg = direction, choices=c("up", "down"))

  purrr::map(seq_along(deg_res), function(i){
    degs <-
      deg_res[[i]] %>%
        tibble::as_tibble(
          rownames = "gene"
        )

    degs <-
      switch(
        EXPR = direction,
        up   =
          dplyr::filter(
            .data = degs,
            !is.na(padj) & padj <= 0.05,
            log2FoldChange > 0
            ),
        down =
          dplyr::filter(
            .data = degs,
            !is.na(padj) & padj <= 0.05,
            log2FoldChange < 0
          )
      )

    degs %>%
      dplyr::mutate(
        log2FoldChange = abs(log2FoldChange),
        dplyr::across(
          .cols = tidyselect::where(is.numeric),
          .fns = \(x) signif(x = x, digits =  2)
          )
        ) %>%
      dplyr::top_n(
        n = 25,
        wt = log2FoldChange
        ) %>%
      dplyr::arrange(
        dplyr::desc(
          log2FoldChange
          )
        )
    }) %>%
    rlang::set_names(
      nm = purrr::map_chr(
        .x = comparison_list,
        .f = str_remove,
        pattern = stringr::str_glue("{grouping_variable}_")
      )
    ) %>%
    purrr::keep(~ nrow(.x) > 0)
}


# create_upregulation_tables <- function(
#   results,
#   comparison_list,
#   grouping_variable
# ){
#   map(seq_along(results), function(i){
#     results[[i]] %>%
#     as_tibble(
#       rownames = "gene"
#     ) %>%
#     filter(
#       !is.na(padj) & padj <= 0.05,
#       log2FoldChange > 0
#     ) %>%
#     mutate_at(
#       vars(-gene),
#       list(~signif(., 2)
#       )
#     ) %>%
#     top_n(
#       n = 25,
#       wt = log2FoldChange
#     ) %>%
#     arrange(
#       desc(
#         log2FoldChange
#       )
#     )
# }) %>%
#   set_names(
#     nm = map_chr(
#       .x = comparison_list,
#       .f = str_remove,
#       pattern = str_glue("{grouping_variable}_")
#     )
#   ) %>%
#   keep(~ nrow(.x) > 0)
# }

extract_de_genes <- function(
  results,
  comparison_list,
  grouping_variable
){
  purrr::map(
    seq_along(results),
    function(i){
      results[[i]] %>%
        tibble::as_tibble(rownames = "gene") %>%
        dplyr::filter(padj < 0.05) %>%
        dplyr::filter(abs(log2FoldChange) >= 0.5) %>%
        dplyr::pull(gene)
      }) %>%
  rlang::set_names(
    nm = comparison_list
  )
}


group_degs <- function(degs, comparison_vars){
  comparison_vars_regex = paste0(comparison_vars, "_", collapse="|")

  tibble::enframe(degs) %>%
    tidyr::unnest(cols = c(value)) %>%
    dplyr::transmute(
      group =
        stringr::str_extract(
          string = name,
          pattern = comparison_vars_regex
          ) %>%
        stringr::str_remove(pattern = "_$"),
      comparison =
        stringr::str_remove(
          string = name,
          pattern = comparison_vars_regex
        ),
      gene_symbol = value
      ) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(
      gene_symbol = gene_symbol,
      count = n(),
      comparison =
        dplyr::case_when(
          count == 1 ~ comparison,
          count > 1 ~ "multiple"
        )
    ) %>%
    dplyr::select(-count) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames(var = "gene_symbol")
}


calc_deg_means <- function(
  exprs,
  deg_class,
  metadata,
  grouping_variable
){
  if(is.character(grouping_variable)){
    diffused_grouping_variable <- sym(grouping_variable)
  }

  diffused_grouping_variable <- enquo(diffused_grouping_variable)

  deg_means <-
    exprs %>%
    t() %>%
    as_tibble(rownames = "name") %>%
    select(name, one_of(deg_class[["gene_symbol"]])) %>%
    pivot_longer(
      -name,
      names_to = "gene",
      values_to = "expr"
    ) %>%
    left_join(
      metadata
    ) %>%
    group_by(
      gene,
      {{diffused_grouping_variable}}
    ) %>%
    summarise(avg = mean(expr)) %>%
    pivot_wider(
      names_from = gene,
      values_from = avg
      ) %>%
    column_to_rownames(grouping_variable)
}


extract_top_degs <- function(up_tables, down_tables){
  up_degs <-
    purrr::map(
      .x = names(up_tables),
      .f = function(i){
        up_tables[[i]] %>%
          dplyr::filter(padj < 0.05) %>%
          dplyr::top_n(25, log2FoldChange) %>%
          dplyr::pull(gene)
        }
      ) %>%
    rlang::set_names(names(up_tables)) %>%
    unlist() # unlisting a map - should we just use `walk` here instead?

  down_degs <-
    purrr::map(
      .x = names(up_tables),
      .f = function(i){
        up_tables[[i]] %>%
          dplyr::filter(padj < 0.05) %>%
          dplyr::top_n(25, log2FoldChange) %>%
          dplyr::pull(gene)
        }
      ) %>%
    rlang::set_names(names(up_tables)) %>%
    unlist()

  unique(c(up_degs, down_degs))
  }
