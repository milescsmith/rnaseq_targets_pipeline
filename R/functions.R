`%nin%` <- purrr::compose(`!`, `%in%`)

deduplicate_samples <- function(md, samples){
  if (nrow(get_dupes(md, sample_name)) > 0){
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples)) %>%
      mutate(sample_name = make_clean_names(string = sample_name,
                                            case = "all_caps"))
    deduplicated_samples = `names<-`(samples, make_clean_names(string = names(samples),
                                                                        case = "all_caps"))
  } else {
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples))
    deduplicated_samples = samples
  }

  return(list(md = deduplicated_md,
              samples = deduplicated_samples))
}


remove_outliers <- function(dds,
                            pc1_zscore_cutoff,
                            pc2_zscore_cutoff){
  dds  <- estimateSizeFactors(dds,
                              locfun = genefilter::shorth,
                              type = "poscounts")
  vsd <- assay(vst(dds))
  pca_res = irlba::prcomp_irlba(vsd)[['rotation']] %>%
    as_tibble() %>%
    mutate(sample = colnames(vsd),
           pc1_zscore = abs((PC1 - mean(PC1))/sd(PC1)),
           pc2_zscore = abs((PC2 - mean(PC2))/sd(PC2))) %>%
    inner_join(as_tibble(colData(dds), rownames = "sample"))

  pc1_outliers <- pca_res %>%
    filter(pc1_zscore >= pc1_zscore_cutoff) %>%
    pull(sample)

  pc2_outliers <- pca_res %>%
    filter(pc2_zscore >= pc2_zscore_cutoff) %>%
    pull(sample)

  outliers <- unique(c(pc1_outliers, pc2_outliers))

  if (length(outliers > 0)){
    dds <- dds[,colnames(dds) %nin% outliers]
  }
  return(list(dds = dds, pca = pca_res, removed = outliers))
}


#' #' optimize_cluster_resolution
#' #'
#' #' @param snn_graph SNN graph output from \code{Seurat::FindNeighbors}
#' #' @param dist_mat Distance matrix from \code{dist()}
#' #' @param from Starting resolution to test
#' #' @param to Final resolution to test
#' #' @param by Amount to increase the resolution each iteration
#' #'
#' #' @return tibble with the silhouette coefficient
#' #' for each resolution tested.
#' #' @export
#' #'
#' #' @examples
#' optimize_cluster_resolution <- function(snn_graph,
#'                                         dist_mat,
#'                                         from = 0.2,
#'                                         to = 1,
#'                                         by = 0.1){
#'
#'
#'   cluster_opt <- future_map_dfc(seq(from = from,
#'                                            to = to,
#'                                            by = by),
#'                                        function(i){
#'                                          leiden(snn_graph,
#'                                                 n_iterations = 10,
#'                                                 resolution_parameter = i,
#'                                                 seed = 8309,
#'                                                 partition_type = "RBConfigurationVertexPartition")
#'                                        })
#'
#'   cluster_opt <- cluster_opt %>%
#'     mutate(sample_name = rownames(snn_graph)) %>%
#'     select(sample_name, everything())
#'
#'   cluster_silhouette <-
#'     future_map_dbl(names(select(cluster_opt, -sample_name)),
#'                           function(i){
#'                             sil <- silhouette(
#'                               deframe(
#'                                 select_at(
#'                                   cluster_opt,
#'                                   vars("sample_name", i))
#'                               ),
#'                               dist_mat
#'                             )
#'
#'                             if(!is.na(sil)){
#'                               return(summary(sil)[["si.summary"]][["Mean"]])
#'                             } else {
#'                               return(0)
#'                             }
#'
#'                           })
#'
#'   cs <- tibble(res = seq(from = from,
#'                          to = to,
#'                          by = by),
#'                coeff = cluster_silhouette)
#'
#'   return(cs)
#' }


###--- kmeans version ---###
sample_clustering <- function(
  exprs_mat,
  from = 2,
  to = 20,
  by = 1){

  kmeans_clusters <-
    future_map(seq(from = from,
                   to = to,
                   by = by),
               function(i){
                 clara(
                   x = exprs_mat,
                   k = i,
                   metric = "jaccard",
                   stand = TRUE,
                   samples = 50,
                   pamLike = TRUE
                   )
                 }
               ) %>%
    set_names(seq(from = from,
                  to = to,
                  by = by))

  sils <- map(kmeans_clusters, function(j){
    silhouette(j)
  })

  avg_sils <- map_dbl(sils, function(k){
    mean(k[,3])
  }) %>%
    set_names(seq(from = from,
                  to = to,
                  by = by)) %>%
    enframe() %>%
    mutate(name = as.integer(name))

  optimal_k =
    avg_sils %>%
    filter(name > 2) %>%
    top_n(
      n = 1,
      wt = value
    )

  clusters <- kmeans_clusters[[as.character(optimal_k[["name"]])]][["clustering"]] %>%
    enframe(name = "sample_name",
            value = "cluster") %>%
    mutate(cluster = as_factor(cluster))

  return(list(avg_sils = avg_sils,
              optimal_k = optimal_k,
              clusters = clusters))
}

# Use the output from cluster_silhouette()
plot_resolution_silhouette_coeff <- function(cluster_silhouette){
  cluster_silhouette %>%
    ggplot(aes(x = res,
               y = coeff)) +
    geom_line() +
    theme_cowplot()
}

plot_dispersion_estimate <- function(object, CV = FALSE){
  px <- mcols(object)$baseMean
  sel <- (px > 0)
  px <- px[sel]
  f <- ifelse(CV, sqrt, I)
  py <- f(mcols(object)$dispGeneEst[sel])
  ymin <- 10^floor(log10(min(py[py > 0], na.rm = TRUE)) -
                     0.1)

  outlier_shape <- ifelse(mcols(object)$dispOutlier[sel],
                          1, 16)
  outlier_size <- ifelse(mcols(object)$dispOutlier[sel],
                         2 * 0.45, 0.45)
  outlier_halo <- ifelse(mcols(object)$dispOutlier[sel],
                         "final", "gene-est")

  disp_data <- tibble(px = px,
                      py = pmax(py, ymin),
                      outlier_shape = as_factor(outlier_shape),
                      outlier_size = as_factor(outlier_size),
                      outlier_halo = as_factor(outlier_halo),
                      dispersions = f(dispersions(object)[sel]),
                      dispersions_fit = f(mcols(object)$dispFit[sel]))

  disp_plot <- disp_data %>%
    ggplot(aes(x = px,
               y = py)) +
    geom_point() +
    geom_point(aes(x = px,
                   y = dispersions,
                   size = outlier_size,
                   shape = outlier_shape,
                   color = outlier_halo)) +
    scale_x_log10() +
    scale_y_log10() +
    scale_shape_manual(values = c(1, 16)) +
    scale_size_manual(values = c(1,2)) +
    scale_color_manual(values = c(
      "dodgerblue",
      "red",
      "black"), ) +
    geom_line(mapping = aes(x = px,
                            y = dispersions_fit,
                            color = "fitted"),
              size = 1) +
    labs(x = "mean of normalized counts",
         y = "dispersion",
         color = "") +
    guides(size = "none",
           shape = "none") +
    theme_cowplot() +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

  return(disp_plot)
}

fix_antibody_values <- function(i) {
  recode(.x = i,
         Negative = "negative",
         POSITIVE = "positive",
         `no_val` = "no_val",
         Indeterminate = "indeterminate")
}


#' Title
#'
#' @param object DESeqResults object
#'
#' @return
#' @export
#'
#' @examples
alt_summary <- function(object){
  notallzero <- sum(object$baseMean > 0)
  up <- sum(object[["padj"]] < 0.05 & object$log2FoldChange >
              metadata(object)$lfcThreshold, na.rm = TRUE)
  down <- sum(object[["padj"]] < 0.05 & object$log2FoldChange <
                metadata(object)$lfcThreshold, na.rm = TRUE)
  outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
  if (is.null(metadata(object)$filterThreshold)) {
    ft <- 0
  } else {
    ft <- round(metadata(object)$filterThreshold)
  }

  filt <- sum(!is.na(object$pvalue) & is.na(object$padj))

  total <- nrow(object)


  tibble(
    up = up,
    down = down,
    outlier = outlier,
    ft = ft,
    lowcounts = filt,
    total = total
  )
}
