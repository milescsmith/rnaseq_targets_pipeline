create_results_list <- function(
  comparison_list,
  dds,
  comparison_grouping_variable
  ){
  map(comparison_list, function(i) {
    lfcShrink(
      dds = dds,
      coef = i,
      parallel = TRUE,
      type = "apeglm"
      )
  }) %>%
    set_names(
      map_chr(
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
  direction <- match.arg(direction, choices=c("up", "down"))

  map(seq_along(deg_res), function(i){
    degs <-
      deg_res[[i]] %>%
        as_tibble(
          rownames = "gene"
        )

    if (direction == "up"){
      degs <-
        filter(
          .data = degs,
          !is.na(padj) & padj <= 0.05,
          log2FoldChange > 0
        )
    } else if (direction == "down"){
      degs <-
        filter(
          .data = degs,
          !is.na(padj) & padj <= 0.05,
          log2FoldChange < 0
        )
    }

    degs %>%
      mutate(log2FoldChange = abs(log2FoldChange)) %>%
      mutate_at(
        .vars = vars(-gene),
        .funs = list(~signif(x = ., digits =  2))
      ) %>%
      top_n(
        n = 25,
        wt = log2FoldChange
      ) %>%
      arrange(
        desc(
          log2FoldChange
        )
      )
  }) %>%
    set_names(
      nm = map_chr(
        .x = comparison_list,
        .f = str_remove,
        pattern = str_glue("{grouping_variable}_")
      )
    ) %>%
    keep(~ nrow(.x) > 0)
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
  map(seq_along(results), function(i){
    results[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        padj < 0.05
      ) %>%
      filter(
        abs(
          log2FoldChange
        ) >= 0.5
      ) %>%
      pull(gene)
  }) %>%
  set_names(
    map_chr(
      .x = comparison_list,
      .f = str_remove,
      pattern = str_glue("{grouping_variable}_")
    )
  )
}


group_degs <- function(degs){
  enframe(degs) %>%
    unnest(cols = c(value)) %>%
    group_by(value) %>%
    mutate(count = n()) %>%
    mutate(
      name =
        case_when(
          count == 1 ~ name,
          count > 1 ~ "multiple"
          )
      ) %>%
    select(-count) %>%
    distinct() %>%
    column_to_rownames(var = "value") %>%
    set_names("comparison")
}


calc_deg_means <- function(
  exprs,
  deg_class,
  metadata,
  grouping_variable
){
  grouping_variable = enquo(grouping_variable)

  deg_means <-
    exprs %>%
    t() %>%
    as_tibble(rownames = "name") %>%
    select(name, one_of(rownames(deg_class))) %>%
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
      {{grouping_variable}}
    ) %>%
    summarise(avg = mean(expr)) %>%
    pivot_wider(
      names_from = gene,
      values_from = avg
      ) %>%
    column_to_rownames("study_group")
}

