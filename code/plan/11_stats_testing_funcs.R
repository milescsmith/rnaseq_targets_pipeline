
modules_compare_with_stats <- function(module_score_table, compare_by){
  module_score_table %>%
    group_by(module) %>%
    wilcox_test(
      as.formula(str_glue("score ~ {compare_by}")),
      p.adjust.method = "BH"
    ) %>%
    grouped_add_xy_positions(
      stats_tbl       = .,
      data_tbl        = module_score_table,
      group_var       = module,
      compare_value   = score
    )
}

pivot_module_scores <- function(module_scores, cluster_var, grouping_var){
  if (is.character(grouping_var)){
    grouping_sym <- sym(grouping_var)
    grouping_sym <- enquo(grouping_sym)
  } else {
    grouping_sym <- enquo(grouping_var)
  }

  if (is.character(cluster_var)){
    cluster_sym <- sym(cluster_var)
    cluster_sym <- enquo(cluster_sym)
  } else {
    cluster_sym <- enquo(cluster_var)
  }

  mutate(
    .data = module_scores,
    {{cluster_sym}}  := as_factor({{cluster_sym}}),
    {{grouping_sym}} := as_factor({{grouping_sym}})
  ) %>%
  select(
    sample_name,
    {{cluster_sym}},
    {{grouping_sym}},
    matches("^M[[:digit:]]+"),
    mg,
    starts_with("ldg")
  ) %>%
  pivot_longer(
    cols = c(
      matches("^M[[:digit:]]+"),
      mg,
      starts_with("ldg")
    ),
    names_to      = "module",
    values_to     = "score"
  )
}
