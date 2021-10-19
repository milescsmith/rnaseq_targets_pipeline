ident_clusters <- function(
  expr_mat,
  sig_pc_method = c("elbow", "horn"),
  max_k = 10,
  ...
){
  sig_pc_method = match.arg(sig_pc_method)

  expr_mat <- expr_mat[which(apply(expr_mat, 1, var) != 0),]

  message("Performing PCA on expression matrix")
  pca_res <-
    PCAtools::pca(
      mat    = expr_mat,
      center = TRUE,
      scale  = FALSE,
      removeVar = 0.1
    )

  pcs_use <-
    switch(
      sig_pc_method,
      elbow = PCAtools::findElbowPoint(pca_res$variance),
      horn = PCAtools::parallelPCA(mat = expr_mat)[["n"]]
    )


  cbs <- future.apply::future_lapply(
    X = seq(2,max_k),
    FUN = \(y) {
      fpc::clusterboot(
        data          = pca_res[["rotated"]][,1:pcs_use],
        clustermethod = claraCBI,
        k             = y
      )
    }
  )

  # TODO: handle situations where there are no
  # identifiable stable clusters
  optimal_k <-
    which(
      purrr::imap_lgl(
        .x = cbs,
        .f = \(x, y) {
          all(cbs[[y]][["bootmean"]] > 0.6)
          }
        )
      ) |>
    
    if (length(optimal_k == 0)){
      sample_clusters <- 
        tibble(
          sample_name = colnames(expr_mat)
          cluster     = 1
        )
    } else {
      sample_clusters <-
        cbs[[max(optimal_k)]][["partition"]] |>
        rlang::set_names(colnames(expr_mat)) |>
        tibble::enframe(
          name  = "sample_name",
          value = "cluster"
        ) %>%
        dplyr::mutate(
          cluster = forcats::as_factor(cluster)
        )
    }

  ret_values <-
    list(
      kmeans_res     = cbs,
      k              = max(optimal_k),
      clusters       = sample_clusters
    )

  ret_values
}
