# Variation can also be examined in reduced dimensional space by PCA or UMAP:
run_pca <- function(
  expr_data,
  metadata
){
  prcomp_irlba(x = expr_data) %>%
    pluck("rotation") %>%
    as_tibble() %>%
    mutate(sample_name = colnames(expr_data)) %>%
    inner_join(metadata)
}

run_umap <- function(
  expr_data,
  metadata
  ){
  umap_results =
    umap(
      t(expr_data),
      n_threads = detectCores(),
      n_sgd_threads = detectCores(),
      verbose = TRUE,
      n_components = 3,
      n_neighbors = 5,

    ) %>%
    as_tibble(.name_repair = "unique") %>%
    set_names(c("umap_1",
                "umap_2",
                "umap_3")) %>%
    mutate(sample_name = colnames(expr_data)) %>%
    inner_join(metadata)
}
