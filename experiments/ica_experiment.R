# TODO: make an optional accessory gene list report
# tar_target(
#   name = plafker_gene_list_file,
#   "references/plafker_gene_list.csv",
#   format = "file"
# ),

# tar_target(
#   name = plafker_gene_list,
#   command = create_module_list(plafker_gene_list_file)
# ),

# tar_target(
#   name = pathway_exprs,
#   command =
#     map(plafker_gene_list, extract_pathway_exprs, exprs_mat = vsc_exprs, metadata = study_md),
# ),

# tar_target(
#   name = pathway_expr_stats,
#   command = map(pathway_exprs, calc_gex_stats)
# ),

# tar_target(
#   name = dds_with_pathway_scores,
#   command =
#     scoreEigengenes(
#       object = dataset_with_scores,
#       module_list = plafker_gene_list,
#       score_func = 'rsvd'
#     )
# ),

# tar_target(
#   name = pathway_eigenvalues,
#   command = extract_module_scores(dds_with_pathway_scores, names(plafker_gene_list))
# ),

# tar_target(
#   name = pathway_eigenvalues_long,
#   command =
#     pathway_eigenvalues %>%
#     pivot_longer(
#       -sample_name,
#       names_to      = "pathway",
#       values_to     = "score"
#     ) %>%
#     left_join(
#       select(
#         .data       = tar_read(study_md),
#         sample_name,
#         study_group
#       )
#     )
#   ),

# tar_target(
#   name = pathway_eigenvalues_stats,
#   command =
#     pathway_eigenvalues_long %>%
#     group_by(pathway) %>%
#     wilcox_test(
#       score ~ study_group,
#       ref.group = "control"
#     ) %>%
#     adjust_pvalue(method = "BH") %>%
#     add_significance() %>%
#     grouped_add_xy_positions(
#       stats_tbl = .,
#       data_tbl = pathway_eigenvalues_long,
#       group_var = pathway,
#       compare_value = score,
#       percent_shift_down = 0.99
#     )
# ),

# This... did not work out so great.  ICA not the best for figuring out the modules
# but this could be adapted for PLIER or NMF or whatever
# tar_target(
#   name = ica_res,
#   command =
#     icafast(
#       X = vsc_top,
#       nc = 25,
#       alg = "par"
#       ),
#   packages = "ica"
# ),
#
# tar_target(
#   name = ica_eigenvalues,
#   command =
#     pluck(
#       .x = ica_res,
#       "S"
#     ) %>%
#     as.data.frame() %>%
#     set_rownames(rownames(vsc_top)) %>%
#     set_colnames(str_glue("IC{1:25}"))
# ),
#
# tar_target(
#   name = ica_eigenvalues_long,
#   command =
#     as_tibble(
#       x = ica_eigenvalues,
#       rownames = "sample_name") %>%
#     pivot_longer(
#       cols = starts_with("IC"),
#       names_to = "IC",
#       values_to = "score"
#       ) %>%
#     left_join(
#       y =
#         select(
#           .data = study_md,
#           sample_name,
#           cluster,
#           study_group
#           )
#       )
# ),
#
# tar_target(
#   name = ica_eigenvalues_stats,
#   command =
#     group_by(
#       .data = ica_eigenvalues_long,
#       pathway
#       ) %>%
#     wilcox_test(
#       score ~ study_group,
#       ref.group = "control"
#       ) %>%
#     adjust_pvalue(method = "BH") %>%
#     add_significance() %>%
#     grouped_add_xy_positions(
#       stats_tbl = .,
#       data_tbl = ica_eigenvalues_long,
#       group_var = IC,
#       compare_value = score,
#       percent_shift_down = 0.99
#       )
#   ),
#
# tar_target(
#   name = ica_loadings,
#   command =
#     pluck(
#       .x = ica_res,
#       "M"
#       ) %>%
#     as.data.frame() %>%
#     set_rownames(colnames(vsc_top)) %>%
#     set_colnames(str_glue("IC{1:25}"))
#   ),
#
# tar_target(
#   name = ic_top_25_genes,
#   command =
#     ica_loadings %>%
#     as_tibble(rownames = "gene") %>%
#     pivot_longer(
#       -gene,
#       names_to = "IC",
#       values_to = "score"
#     ) %>%
#     group_by(IC) %>%
#     top_n(
#       n = 25,
#       wt = score
#     ) %>%
#     select(-score) %>%
#     rename(module = "IC") %>%
#     mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
#     filter(!is.na(hugo))
#   ),
#
# tar_target(
#   name = ic_enriched_list,
#   command =
#     map_dfr(
#       .x = unique(ic_top_25_genes$module),
#       .f = function(x){
#         module_gsea(
#           module_genes = ic_top_25_genes,
#           module_of_interest = x)
#       }
#     ),
#   packages = "purrr"
#   ),
#
# tar_target(
#   name = ic_plotting,
#   command =
#     module_gsea_plots(ic_enriched_list) %>%
#     mutate(
#       module =
#         str_remove(
#           string = module,
#           pattern = "^ME"
#           )
#       )
# ),
