# Variation can also be examined in reduced dimensional space by PCA or UMAP:
run_pca <- function(
  expr_data,
  metadata
){
  irlba::prcomp_irlba(x = expr_data) |>
    purrr::pluck("rotation") |>
    tibble::as_tibble() |>
    dplyr::mutate(sample_name = colnames(expr_data)) |>
    dplyr::inner_join(metadata)
}

run_umap <- function(
  expr_data,
  metadata
  ){
  uwot::umap(
    t(expr_data),
    n_threads = parallel::detectCores(),
    n_sgd_threads = parallel::detectCores(),
    verbose = TRUE,
    n_components = 3,
    n_neighbors = 5
  ) |>
  tibble::as_tibble(.name_repair = "unique") |>
  rlang::set_names(
    c(
      "umap_1",
      "umap_2",
      "umap_3"
      )
    ) |>
  dplyr::mutate(sample_name = colnames(expr_data)) |>
  dplyr::inner_join(metadata)
}
