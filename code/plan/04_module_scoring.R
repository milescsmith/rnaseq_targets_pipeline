
dds_with_scores =
  scoreEigengenes(
    object = dds_processed,
    module_list = banchereau_modules,
    score_func = 'rsvd'
    ) %>%
  scoreEigengenes(
    object = .,
    module_list = ldg_modules,
    score_func = 'rsvd'
    ) %>%
  scoreEigengenes(
    object = .,
    module_list = metasignature_module,
    score_func = 'rsvd'
    )

#### The tirosh_score_modules() func isn't working at the moment - something
#### about not enough genes to cut into 100 groups?
# banchereau_module_scores =
#   tirosh_score_modules(
#     expr_obj = vsd_exprs,
#     module_list = banchereau_modules
#   ) %>%
#   as_tibble(rownames = "sample_name")
#
# ldg_module_scores =
#   tirosh_score_modules(
#     expr_obj = vsd_exprs,
#     module_list = ldg_modules
#   ) %>%
#   as_tibble(rownames = "sample_name")
#
# metasig_scores =
#   tirosh_score_modules(
#     vsd_exprs,
#     metasignature_module
#   ) %>%
#   as_tibble(rownames = "sample_name")

# module_scores =
#   inner_join(
#     banchereau_module_scores,
#     ldg_module_scores
#   ) %>%
#   inner_join(metasig_scores) %>%
#   column_to_rownames(var = "sample_name")

module_scores =
  colData(dds_with_scores) %>%
  as_tibble(rownames = "sample_name") %>%
  select(
    sample_name,
    one_of(
      c(
        names(banchereau_modules),
        names(ldg_modules),
        names(metasignature_module)
        )
      )
    )

# dds_with_scores = target({
#   colData(dds_processed) <-
#     left_join(
#       as_tibble(
#         colData(dds_processed),
#         rownames = "sample"
#       ),
#       as_tibble(
#         module_scores,
#         rownames = "sample"
#       )
#     ) %>%
#     column_to_rownames(var = "sample") %>%
#     DataFrame()
#
#   dds_processed
# })

disp_plot = plot_dispersion_estimate(dds_with_scores)

annotated_modules =
  module_annotation %>%
  filter(type != "Undetermined")

annotated_mod_list =
  annotated_modules %>%
  mutate(
    module_type =
      paste(
        module,
        type,
        sep = " - "
      )
  ) %>%
  select(-type) %>%
  deframe()

annotated_module_scores =
  module_scores %>%
  select(sample_name, one_of(annotated_modules$module))

module_tbl =
  c(
    banchereau_modules,
    ldg_modules,
    metasignature_module
  ) %>%
  enframe() %>%
  unnest(cols = "value") %>%
  rename(
    module = name,
    gene = value
    ) %>%
  inner_join(module_annotation)

module_ISGs = module_tbl %>%
  dplyr::group_by(gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::filter(
    type == "Interferon",
    gene %in% rownames(vsd_exprs)
    ) %>%
  dplyr::select(module, gene) %>%
  arrange(module) %>%
  column_to_rownames("gene")
