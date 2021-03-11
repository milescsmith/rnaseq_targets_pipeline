comparison_results_list =
  resultsNames(dds_with_scores) %>%
  keep(
    str_detect(
      string = .,
      pattern = comparison_grouping_variable
    )
  )

res = map(comparison_results_list, function(i) {
  lfcShrink(dds_with_scores,
            coef = i,
            parallel = TRUE,
            type = "apeglm")
}) %>%
  set_names(
    map_chr(
      comparison_results_list,
      str_remove,
      pattern = paste0(comparison_grouping_variable, "_")
    )
  )

down_tables = map(seq_along(res), function(i){
  res[[i]] %>%
    as_tibble(
      rownames = "gene"
    ) %>%
    filter(
      !is.na(padj) & padj <= 0.05,
      log2FoldChange < 0
    ) %>%
    mutate(
      log2FoldChange = -log2FoldChange
    ) %>%
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
      .x = comparison_results_list,
      .f = str_remove,
      pattern = str_glue("{comparison_grouping_variable}_")
    )
  ) %>%
  keep(~ nrow(.x) > 0)

up_tables = map(seq_along(res), function(i){
  res[[i]] %>%
    as_tibble(
      rownames = "gene"
    ) %>%
    filter(
      !is.na(padj) & padj <= 0.05,
      log2FoldChange > 0
    ) %>%
    mutate_at(
      vars(-gene),
      list(~signif(., 2)
      )
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
      .x = comparison_results_list,
      .f = str_remove,
      pattern = str_glue("{comparison_grouping_variable}_")
    )
  ) %>%
  keep(~ nrow(.x) > 0)

degs = map(seq_along(res), function(i){
  res[[i]] %>%
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
      .x = comparison_results_list,
      .f = str_remove,
      pattern = str_glue("{comparison_grouping_variable}_")
    )
  )

deg_class = enframe(degs) %>%
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

deg_means =
  vsd_exprs %>%
  t() %>%
  as_tibble(rownames = "name") %>%
  select(name, one_of(rownames(deg_class))) %>%
  pivot_longer(-name,
               names_to = "gene",
               values_to = "expr") %>%
  left_join(
    as_tibble(
      colData(dds_with_scores),
      rownames = "name"
    ) %>%
      select(
        name,
        sex,
        ethnicity,
        disease_class
      )
  ) %>%
  group_by(gene,
           disease_class) %>%
  summarise(avg = mean(expr)) %>%
  pivot_wider(names_from = gene,
              values_from = avg) %>%
  column_to_rownames("disease_class")

ISGs = intersect(
  c(
    "STAT1", "ADAR", "ABCE1", "RNASEL", "TYK2", "IFNAR1",
    "IFNB1", "STAT2", "IFNAR2", "JAK1", "SAMHD1",
    "SOCS1", "SOCS3", "STAT1", "ISG20", "IFITM3",
    "IFITM1", "IRF9", "ISG15", "IFI6", "IFIT3", "USP18",
    "IP6K2", "PSMB8", "IFIT1", "IRF4", "IRF5", "IRF1", "IRF3",
    "IRF6", "IRF8", "IRF2", "IRF7", "IFITM2", "XAF1", "IFI27",
    "GBP2", "RSAD2", "MX2", "MX1", "IFIT2", "IFI35", "BST2",
    "OAS1", "OASL", "OAS2", "OAS3", "PTPN6", "PTPN11"
    ),
  rownames(vsd_exprs)
  )
