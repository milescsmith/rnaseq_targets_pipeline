`%nin%` <- compose(`!`, `%in%`)

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


remove_outliers <- function(dds, zscore_cutoff){
  dds  <- estimateSizeFactors(dds,
                              locfun = shorth,
                              type = "poscounts")
  vsd <- assay(vst(dds))
  pca_res = prcomp_irlba(vsd)[['rotation']] %>%
    as_tibble() %>%
    mutate(sample = colnames(vsd),
           zscore = abs((PC1 - mean(PC1))/sd(PC1))) %>%
    inner_join(as_tibble(colData(dds), rownames = "sample"))
  
  outliers <- pca_res %>% filter(zscore >= zscore_cutoff) %>% pull(sample)
  if (length(outliers > 0)){
    dds <- dds[,colnames(dds) %nin% outliers] 
  }
  return(list(dds = dds, pca = pca_res, removed = outliers))
}


#' optimize_cluster_resolution
#'
#' @param snn_graph SNN graph output from \code{Seurat::FindNeighbors}
#' @param dist_mat Distance matrix from \code{dist()}
#' @param from Starting resolution to test
#' @param to Final resolution to test
#' @param by Amount to increase the resolution each iteration
#'
#' @return tibble with the silhouette coefficient
#' for each resolution tested.
#' @export
#'
#' @examples
optimize_cluster_resolution <- function(snn_graph,
                                        dist_mat,
                                        from = 0.2,
                                        to = 1,
                                        by = 0.1){
  
  cluster_opt <- furrr::future_map_dfc(seq(from = from,
                                           to = to,
                                           by = by),
                                       function(i){
                                         leiden(as.matrix(snn_graph),
                                                n_iterations = 10,
                                                resolution_parameter = i,
                                                seed=as.numeric(Sys.time()),
                                                partition_type = "RBConfigurationVertexPartition")
                                       })
  
  cluster_opt <- cluster_opt %>%
    mutate(sample_name = rownames(snn_graph)) %>%
    select(sample_name, everything())
  
  cluster_silhouette <- 
    furrr::future_map_dbl(names(select(cluster_opt, -sample_name)),
                          function(i){
                            sil <- silhouette(
                              deframe(
                                select_at(
                                  cluster_opt,
                                  vars("sample_name", i))
                              ),
                              dist_mat
                            )
                            
                            if(!is.na(sil)){
                              return(summary(sil)[["si.summary"]][["Mean"]])
                            } else {
                              return(0)
                            }
                            
                          })
  
  cs <- tibble(res = seq(from = from,
                         to = to,
                         by = by),
               coeff = cluster_silhouette)
  
  return(cs)
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