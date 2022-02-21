save_vsd_exprs =
  vsd_exprs %>%
  as.data.frame() %>%
  fwrite(
    row.names = TRUE,
    file = file_out("results/filtered_normalized_stabilized_expression.csv")
  )

save_study_md =
  study_md %>%
  fwrite(
    row.names = TRUE,
    file = file_out("results/filtered_metadata.csv")
  )

save_counts =
  counts(dds_with_scores) %>%
  as.data.frame() %>%
  fwrite(
    row.names = TRUE,
    file = file_out("results/filtered_transcript_counts.csv")
  )

save_results =
  res %>%
  as.data.frame() %>%
  fwrite(
    row.names = TRUE,
    file = file_out("results/log_fold_changes.csv")
  )

save_wgcna_modules =
  wgcna_modules$MEs %>%
  fwrite(
    row.names = TRUE,
    file = file_out("results/wgcna_eigengene_scores.csv")
  )

save_wgcna_genes =
  wgcna_module_genes %>%
  write_csv(path = file_out("results/wgcna_genes.csv"))

save_module_scores =
  study_md %>%
  select(
    sample_name,
    matches("^M[[:digit:]]+"),
    one_of(names(ldg_modules))
  ) %>%
  as.data.frame() %>%
  fwrite(
    row.names = TRUE,
    file = file_out("results/module_scores.csv")
  )

save_objects_for_vis = target({
  save(
    annotation_info,
    final_md,
    group_pal,
    module_scores_with_viral,
    ISGs,
    pca_results,
    degs,
    umap_results,
    vsd_exprs,
    wgcna_modules,
    file = file_out("results/for_shiny_vis.RData"))
    })

