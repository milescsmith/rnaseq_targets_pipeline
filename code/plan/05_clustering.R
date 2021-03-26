sample_dists =
  vsd_exprs %>%
  t() %>%
  parallelDist::parallelDist()

#fig.width=12, fig.height=9
sampleDistMatrix = as.matrix(sample_dists)

sample_dendrogram =
  sample_dists %>%
  hclust() %>%
  as.dendrogram()

sample_cluster_info =
  ident_clusters(
    column_to_rownames(annotated_module_scores, "sample_name"),
    K.max = 20
  )

clusters =
  sample_cluster_info$clusters %>% mutate(cluster = as_factor(cluster))

study_md =
  colData(dds_with_scores) %>%
  as_tibble(rownames = "sample_name") %>%
  left_join(clusters)

annotation_info =
  select(.data = study_md,
         disease_class,
         sex,
         cluster,
         sample_name) %>%
  column_to_rownames(var = "sample_name")
