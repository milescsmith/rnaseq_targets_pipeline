#' @title getMetaData
#'
#' @keywords internal
getMetaData <- function(object, ...){
  UseMethod("getMetaData")
}

#' @rdname getMetaData
#' @method getMetaData DESeqDataSet
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble as_tibble
#'
#' @return
#' @keywords internal
getMetaData.DESeqDataSet <-
  function(object, ...){
    SummarizedExperiment::colData(object) |>
      tibble::as_tibble(rownames = "sample_name")
  }

#' @rdname getMetaData
#' @method getMetaData DGEList
#'
#' @importFrom tibble as_tibble
#'
#' @return
#' @keywords internal
getMetaData.DGEList <-
  function(object, ...){
    object[["samples"]] |>
      tibble::as_tibble(rownames = "sample_name")
  }

plot_dispersion_estimate <- function(object,...){
  UseMethod("plot_dispersion_estimate", object)
}

plot_dispersion_estimate.DGEList <- function(object, ...){ NULL }

plot_dispersion_estimate.DESeqDataSet <- function(object, CV = FALSE){
  px <- mcols(object)$baseMean
  sel <- (px > 0)
  px <- px[sel]
  f <- ifelse(CV, sqrt, I)
  py <- f(mcols(object)$dispGeneEst[sel])
  ymin <- 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 0.1)

  outlier_shape <-
    ifelse(
      test = mcols(object)$dispOutlier[sel],
      yes = 1,
      no = 16
    )

  outlier_size <-
    ifelse(
      test = mcols(object)$dispOutlier[sel],
      yes = 2 * 0.45,
      no =  0.45
    )

  outlier_halo <-
    ifelse(
      test = mcols(object)$dispOutlier[sel],
      yes = "final",
      no = "gene-est"
    )

  disp_data <-
    tibble(
      px = px,
      py = pmax(py, ymin),
      outlier_shape = as_factor(outlier_shape),
      outlier_size = as_factor(outlier_size),
      outlier_halo = as_factor(outlier_halo),
      dispersions = f(dispersions(object)[sel]),
      dispersions_fit = f(mcols(object)$dispFit[sel])
    )

  disp_plot <- disp_data %>%
    ggplot(
      aes(
        x = px,
        y = py
      )
    ) +
    geom_point() +
    geom_point(
      aes(
        x = px,
        y = dispersions,
        size = outlier_size,
        shape = outlier_shape,
        color = outlier_halo
      )
    ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_shape_manual(values = c(1, 16)) +
    scale_size_manual(values = c(1,2)) +
    scale_color_manual(values = c(
      "dodgerblue",
      "red",
      "black"), ) +
    geom_line(
      mapping =
        aes(
          x = px,
          y = dispersions_fit,
          color = "fitted"
        ),
      size = 1
    ) +
    labs(
      x = "mean of normalized counts",
      y = "dispersion",
      color = ""
    ) +
    guides(
      size = "none",
      shape = "none"
    ) +
    theme_pubr() +
    theme(
      legend.justification=c(1,0),
      legend.position=c(1,0)
    )

  return(disp_plot)
}


ident_clusters <- function(
  expr_mat,
  optimal_k_method = "Tibs2001SEmax",
  nstart = 25,
  K.max = 50,
  B = 100,
  d.power = 2
){
  module_rf <-
    randomForest::randomForest(
      x = expr_mat,
      y = NULL,
      prox = TRUE)

  rf_distance_mat <-
    stats::dist(1 - module_rf$proximity) %>%
    as.matrix()

  kmeans_gap_stat <-
    cluster::clusGap(
      x          = rf_distance_mat,
      FUNcluster = kmeans,
      nstart     = nstart,
      K.max      = K.max,
      B          = B,
      d.power    = d.power
    )

  new_optimal_k <-
    with(
      data = kmeans_gap_stat,
      expr = cluster::maxSE(
        Tab[,"gap"],
        Tab[,"SE.sim"],
        method=optimal_k_method
      )
    )

  k_clusters <-
    stats::kmeans(
      x       = rf_distance_mat,
      centers = new_optimal_k,
      nstart  = 25
    )

  sample_clusters <-
    tibble::enframe(
      x     = k_clusters[["cluster"]],
      name  = "sample_name",
      value = "cluster"
    ) %>%
    dplyr::mutate(
      cluster = forcats::as_factor(cluster)
    )

  ret_values <-
    list(
      kmeans_res  = k_clusters,
      rf_distance = rf_distance_mat,
      clusters    = sample_clusters,
      gap_stat    = kmeans_gap_stat
    )

  ret_values
}

targets_recode <- function(
  target_list,
  thing_to_unquote_splice
){
  dplyr::recode(
    .x = target_list,
    !!! {{thing_to_unquote_splice}}
  )
}

symbol_or_string <- function(var) {
  v <- enquo(var)
  str(list(
    v = v,
    quo_is_missing = quo_is_missing(v),
    quo_is_null = quo_is_null(v),
    quo_is_symbol = quo_is_symbol(v),
    quo_is_symbolic = quo_is_symbolic(v),
    quo_is_call = quo_is_lang(v)))
}
