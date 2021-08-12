
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

pivot_module_scores <- function(module_scores){
  mutate(
    .data = module_scores,
    responder = as_factor(responder)
  ) %>%
  select(
    sample_name,
    responder,
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
