
pca <- rsvd::rpca(A = t(vsd_exprs), k = 20)
pca_sdev_drop <- c(diff(pca$sdev), 0) / -pca$sdev
sig_pca_dims <- rev(which(pca_sdev_drop > 0.025))[1]

pca_results = irlba::prcomp_irlba(x = vsd_exprs, n = sig_pca_dims) %>%
  `[[`("rotation") %>%
  as_tibble() %>%
  mutate(sample_name = colnames(vsd_exprs)) %>%
  inner_join(as_tibble(colData(dds_processed),
                       rownames="sample_name")) %>%
  inner_join(clusters)


umap_results = uwot::umap(pca_results[,1:6],
                          n_threads = detectCores(),
                          n_sgd_threads = detectCores(),
                          verbose = TRUE,
                          n_components = 3) %>%
  as_tibble(.name_repair = "unique") %>%
  `colnames<-`(c("umap_1",
                 "umap_2",
                 "umap_3")) %>%
  mutate(sample_name = colnames(vsd_exprs)) %>%
  inner_join(as_tibble(colData(dds_processed),
                       rownames="sample_name")) %>%
  inner_join(clusters)

umap_results %>%
  plot_ly(x = ~umap_1,
          y = ~umap_2,
          color = ~cluster,
          alpha = 0.75,
          colors = group_pal$cluster,
          text = ~paste('<br>Sample: ', sample_name,
                        '<br>Classification: ', disease_class,
                        '<br>Ethnicity ', race_code,
                        '<br>Cluster: ', cluster,
                        '<br>Sex: ', sex,
                        '<br>Plate: ', run_id
          )
  ) %>%
  add_markers(marker = list(size = 7.5,
                            line = list(color = 'rgba(0, 0, 0, .8)',
                                        width = 1)))

enframe(kmeans(x = column_to_rownames(umap_results[,1:4],
                                      var = "sample_name"),
               centers = 7)[["cluster"]],
        name = "sample_name",
        value = "k_means") %>%
  inner_join(umap_results) %>%
  mutate(k_means = as_factor(k_means)) %>%
  plot_ly(x = ~umap_1,
          y = ~umap_2,
          z = ~umap_3,
          color = ~k_means,
          size = 2.5,
          alpha = 0.9,
          #colors = group_pal$cluster,
          text = ~paste('<br>Sample: ', sample_name,
                        '<br>Classification: ', disease_class,
                        '<br>Ethnicity ', race_code,
                        '<br>Cluster: ', k_means,
                        '<br>Sex: ', sex,
                        '<br>Plate: ', run_id
          )
  ) %>%
  add_markers(marker = list(size = 3,
                            line = list(color = 'rgba(0, 0, 0, .8)',
                                        width = 1.5)))

enframe(kmeans(x = column_to_rownames(umap_results[,1:4],
                                      var = "sample_name"),
               centers = 5)[["cluster"]],
        name = "sample_name",
        value = "k_means") %>%
  inner_join(umap_results) %>%
  mutate(k_means = as_factor(k_means)) %>%
  plot_ly(x = ~umap_1,
          y = ~umap_2,
          color = ~k_means,
          alpha = 0.75,
          text = ~paste('<br>Sample: ', sample_name,
                        '<br>Classification: ', disease_class,
                        '<br>Ethnicity ', race_code,
                        '<br>Cluster: ', k_means,
                        '<br>Sex: ', sex,
                        '<br>Plate: ', run_id
          )
  ) %>%
  add_markers(marker = list(size = 7.5,
                            line = list(color = 'rgba(0, 0, 0, .8)',
                                        width = 1)))


cluster_opt <- future_map_dfc(seq(from = 1,
                                  to = 50,
                                  by = 1),
                              function(i){
                                kmeans(x = column_to_rownames(umap_results[,1:4],
                                                              var = "sample_name"),
                                       centers = i)[["cluster"]]
                              }) %>% `colnames<-`(paste0("k",1:50))

cluster_opt <- cluster_opt %>%
  mutate(sample_name = pca_results[['sample_name']]) %>%
  select(sample_name, everything())

cluster_silhouette <-
  future_map_dbl(names(select(cluster_opt, -sample_name)),
                 function(i){
                   sil <- silhouette(
                     deframe(
                       select_at(
                         cluster_opt,
                         vars("sample_name", i))
                     ),
                     umap_dist
                   )

                   if(!is.na(sil)){
                     return(summary(sil)[["si.summary"]][["Mean"]])
                   } else {
                     return(0)
                   }

                 })

cs <- tibble(res = seq(from = 1,
                       to = 50,
                       by = 1),
             coeff = cluster_silhouette)

plot_resolution_silhouette_coeff(cs)


silhouette_score <- function(k){
  km <- kmeans(pca_results[,1:18], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(pca_results[,1:18]))
  mean(ss[, 3])
}

k <- 2:20
avg_sil <- future_map_dbl(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

tibble(clusters = k, score = avg_sil) %>%
  ggplot(aes(x = clusters, y = score)) +
  geom_line() +
  theme_cowplot()

