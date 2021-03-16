
modules_compare_with_stats <- function(module_score_table, comparison){
  module_score_table %>%
    group_by(module) %>%
    wilcox_test(
      score ~ {{comparison}},
      p.adjust.method = "BH"
    ) %>%
    grouped_add_xy_positions(
      stats_tbl       = .,
      data_tbl        = annotated_module_scores_pivot,
      group_var       = module,
      compare_value   = score
    )
}

pivot_module_scores <- function(module_scores){
  mutate(
    .data = module_scores,
    cluster       = as_factor(cluster),
    disease_class = as_factor(disease_class)
  ) %>%
  select(
    sample_name,
    cluster,
    disease_class,
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