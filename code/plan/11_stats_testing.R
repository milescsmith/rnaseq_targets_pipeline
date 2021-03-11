annotated_module_scores_with_cluster_class =
  module_scores_with_viral %>%
  mutate(
    cluster = as_factor(cluster)
    # cell_type = as_factor(cell_type)
  ) %>%
  select(
    cluster,
    # cell_type,
    disease_class,
    one_of(annotated_modules$module)
  )

renamed_annotated_module_scores =
  names(annotated_module_scores_with_cluster_class) %>%
  drake_recode(thing_to_unquote_splice = annotated_mod_list) %>%
  set_names(
    nm = .,
    x  = annotated_module_scores_with_cluster_class
  )

annotated_module_scores_pivot =
  renamed_annotated_module_scores %>%
  pivot_longer(
    cols      = starts_with("M"),
    names_to  = "module",
    values_to = "score"
  )

annotated_module_stats_by_cluster =
  annotated_module_scores_pivot %>%
  group_by(module) %>%
  wilcox_test(
    score ~ cluster,
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl       = .,
    data_tbl        = annotated_module_scores_pivot,
    group_var       = module,
    compare_value   = score
  )

# annotated_module_stats_by_cell_type =
#   annotated_module_scores_pivot %>%
#   group_by(module) %>%
#   wilcox_test(
#     score ~ cell_type,
#     p.adjust.method = "BH"
#   ) %>%
#   grouped_add_xy_positions(
#     stats_tbl = .,
#     data_tbl = annotated_module_scores_pivot,
#     group_var = module,
#     compare_value = score
#   )

annotated_module_stats_by_disease =
  annotated_module_scores_pivot %>%
  group_by(module) %>%
  wilcox_test(
    score ~ disease_class,
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl       = .,
    data_tbl        = annotated_module_scores_pivot,
    group_var       = module,
    compare_value   = score
  )

module_scores_pivot =
  module_scores_with_viral %>%
  mutate(
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

module_stats_by_cluster =
  module_scores_pivot %>%
  group_by(module) %>%
  wilcox_test(
    score ~ cluster,
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl       = .,
    data_tbl        = module_scores_pivot,
    group_var       = module,
    compare_value   = score
  )

module_stats_by_disease =
  module_scores_pivot %>%
  group_by(module) %>%
  wilcox_test(
    score ~ disease_class,
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl       = .,
    data_tbl        = module_scores_pivot,
    group_var       = module,
    compare_value   = score
  )

module_scores_with_viral_by_cluster =
  module_scores_with_viral %>%
  select(
    cluster,
    matches("^ME")
    ) %>%
  pivot_longer(
    -cluster,
    names_to     = "module",
    values_to    = "score"
    ) %>%
  mutate(cluster = as_factor(cluster)
  )

module_scores_with_viral_by_cluster_stats =
  module_scores_with_viral_by_cluster %>%
  group_by(module) %>%
  wilcox_test(score ~ cluster) %>%
  grouped_add_xy_positions(
    stats_tbl     = .,
    data_tbl      = module_scores_with_viral_by_cluster,
    group_var     = module,
    compare_value = score
  )

# module_scores_with_viral_by_cell_type =
#   module_scores_with_viral %>%
#   select(cell_type,
#          matches("^ME")) %>%
#   pivot_longer(-cell_type,
#                names_to = "module",
#                values_to = "score") %>%
#   mutate(cell_type = as_factor(cell_type)
#   )

# module_scores_with_viral_by_cell_type_stats =
#   module_scores_with_viral_by_cell_type %>%
#   group_by(module) %>%
#   wilcox_test(score ~ cell_type) %>%
#   grouped_add_xy_positions(
#     stats_tbl = .,
#     data_tbl = module_scores_with_viral_by_cell_type,
#     group_var = module,
#     compare_value = score
#   )

module_scores_with_viral_by_disease =
  module_scores_with_viral %>%
  select(
    disease_class,
    matches("^ME")
    ) %>%
  pivot_longer(
    -disease_class,
    names_to           = "module",
    values_to          = "score"
    ) %>%
  mutate(disease_class = as_factor(disease_class)
  )

module_scores_with_viral_by_disease_stats =
  module_scores_with_viral_by_disease %>%
  group_by(module) %>%
  wilcox_test(score ~ disease_class) %>%
  grouped_add_xy_positions(
    stats_tbl = .,
    data_tbl = module_scores_with_viral_by_disease,
    group_var = module,
    compare_value = score
  )
