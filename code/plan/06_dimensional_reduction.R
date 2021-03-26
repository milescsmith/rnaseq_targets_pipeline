# Variation can also be examined in reduced dimensional space by PCA or UMAP:
pca_results = prcomp_irlba(x = vsd_exprs) %>%
  pluck("rotation") %>%
  as_tibble() %>%
  mutate(sample_name = colnames(vsd_exprs)) %>%
  inner_join(as_tibble(colData(dds_with_scores),
                       rownames="sample_name")) %>%
  inner_join(clusters)

umap_results =
  umap(
    t(vsd_exprs),
    n_threads = detectCores(),
    n_sgd_threads = detectCores(),
    verbose = TRUE,
    n_components = 3
  ) %>%
  as_tibble(.name_repair = "unique") %>%
  set_names(c("umap_1",
              "umap_2",
              "umap_3")) %>%
  mutate(sample_name = colnames(vsd_exprs)) %>%
  inner_join(
    as_tibble(
      colData(dds_with_scores),
      rownames = "sample_name"
    )
  ) %>%
  inner_join(clusters)
