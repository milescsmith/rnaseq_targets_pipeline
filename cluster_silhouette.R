future::plan("multiprocess")

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
                              sample_dists
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

plot_resolution_silhouette_coeff <- function(cluster_silhouette){
  cluster_silhouette %>%
    ggplot(aes(x = res,
               y = coeff)) +
    geom_line() +
    theme_cowplot()
  
}
