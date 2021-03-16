library(caret)
library(yardstick)
library(rsample)
library(tidyverse)

library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 4)

drake::loadd(module_scores, annotation_info)

module_ml_scores <- module_scores %>%
  as_tibble(rownames="sample_name") %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name")) %>%
  select(cluster, starts_with("M"))

module_split <- initial_split(data = module_ml_scores,
                              prop = 0.8,
                              strata = "cluster")

module_train <- training(module_split)
module_test <- testing(module_split)

module_rf_cv <- train(cluster ~ .,
                        method = "parRF",
                        data = module_train,
                        trControl = trainControl(method = "repeatedcv",
                                                 number = 10,
                                                 repeats = 10))

module_rf_boot <- train(cluster ~ .,
                   method = "parRF",
                   data = module_train,
                   trControl = trainControl(method = "boot",
                                            sampling = "up"))

caret::confusionMatrix(
  data = predict(module_rf_cv, module_test),
  reference = module_test$cluster
)

caret::confusionMatrix(
  data = predict(module_rf_boot, module_test),
  reference = module_test$cluster
)

module_rf_cv_varImp <- varImp(object = module_rf_cv, scale = FALSE)

top_classifying_modules <- module_rf_cv_varImp$importance %>%
  as_tibble(rownames='module') %>%
  top_n(10, Overall) %>%
  pull(module)

module_ml_scores %>%
  pivot_longer(-cluster,
               names_to = "module",
               values_to = "score") %>%
  filter(module %in% top_classifying_modules) %>%
  ggplot(aes(x = cluster, y = score)) +
  geom_boxplot() +
  facet_wrap(facets = vars(module), scales = "free_y") +
  theme_cowplot()

#----------------------WGCNA------------------------------#
library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 4)

wgcna_scores <- wgcna_modules$MEs %>%
  as_tibble(rownames="sample_name") %>%
  select(-MEgrey) %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name")) %>%
  select(cluster, starts_with("ME"))

wgcna_split <- initial_split(data = wgcna_scores,
                              prop = 0.8,
                              strata = "cluster")

wgcna_train <- training(wgcna_split)
wgcna_test <- testing(wgcna_split)

wgcna_rf_cv <- train(cluster ~ .,
                      method = "parRF",
                      data = wgcna_train,
                      trControl = trainControl(method = "repeatedcv",
                                               number = 10,
                                               repeats = 10))

wgcna_rf_boot <- train(cluster ~ .,
                        method = "parRF",
                        data = wgcna_train,
                        trControl = trainControl(method = "boot",
                                                 sampling = "up"))

caret::confusionMatrix(
  data = predict(wgcna_rf_cv, wgcna_test),
  reference = wgcna_test$cluster
)

caret::confusionMatrix(
  data = predict(wgcna_rf_boot, wgcna_test),
  reference = wgcna_test$cluster
)

wgcna_rf_cv_varImp <- varImp(object = wgcna_rf_cv,
                             scale = TRUE,
                             useModel = TRUE)

wgcna_scores %>%
  pivot_longer(-cluster,
               names_to = "module",
               values_to = "score") %>%
  ggplot(aes(x = cluster,
             y = score)) +
  geom_boxplot() +
  facet_wrap(facets = vars(module),
             scales = "free_y") +
  theme_cowplot()

vsd_top %>%
  as_tibble(rownames = "sample_name") %>%
  pivot_longer(-sample_name,
               names_to = "gene",
               values_to = "score") %>%
  left_join(as_tibble(annotation_info, rownames = "sample_name")) %>%
  filter(gene %in% chooseTopHubInEachModule(vsd_top, wgcna_modules$colors, power = 4, type = "signed hybrid")) %>%
  ggplot(aes(x = cluster,
             y = score)) +
  geom_boxplot() +
  facet_wrap(facets = vars(gene),
             scales = "free_y") +
  theme_cowplot()
#----------------------------ICA------------------------------#

library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 4)

ica_ml_scores <- ica_scores %>%
  as_tibble(rownames="sample_name") %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name")) %>%
  select(cluster, starts_with("IC"))

ica_split <- initial_split(data = ica_ml_scores,
                             prop = 0.8,
                             strata = "cluster")

ica_train <- training(ica_split)
ica_test <- testing(ica_split)

ica_rf_cv <- train(cluster ~ .,
                     method = "parRF",
                     data = ica_train,
                     trControl = trainControl(method = "repeatedcv",
                                              number = 10,
                                              repeats = 10))

ica_rf_boot <- train(cluster ~ .,
                       method = "parRF",
                       data = ica_train,
                       trControl = trainControl(method = "boot",
                                                sampling = "up"))

ica_xgbDART_cv <- train(cluster ~ .,
                   method = "xgbDART",
                   data = ica_train,
                   trControl = trainControl(method = "repeatedcv",
                                            number = 10,
                                            repeats = 10))

caret::confusionMatrix(
  data = predict(ica_rf_cv, ica_test),
  reference = ica_test$cluster
)

caret::confusionMatrix(
  data = predict(ica_rf_boot, ica_test),
  reference = ica_test$cluster
)

ica_rf_cv_varImp <- varImp(object = ica_rf_cv, scale = FALSE)

top_classifying_ica <- ica_rf_cv_varImp$importance %>%
  as_tibble(rownames='module') %>%
  top_n(10, Overall) %>%
  pull(module)

ica_ml_scores %>%
  pivot_longer(-cluster,
               names_to = "module",
               values_to = "score") %>%
  filter(module %in% top_classifying_ica) %>%
  ggplot(aes(x = cluster, y = score)) +
  geom_boxplot() +
  facet_wrap(facets = vars(module), scales = "free_y") +
  theme_cowplot()

wgcna_scores %>%
  # as_tibble(rownames = "sample_name") %>%
  # inner_join(as_tibble(annotation_info,
  #                      rownames = "sample_name")) %>%
  # inner_join(as_tibble(ica_scores,
  #                      rownames = "sample_name")) %>%
  ggplot(aes(x = MEyellow, y = MEturquoise, color = cluster)) +
  geom_point() +
  cowplot::theme_cowplot()


featurePlot(x = ica_ml_scores[, 2:10],
           y = ica_ml_scores$cluster, 
           plot = "box", 
           ## Pass in options to bwplot() 
           scales = list(y = list(relation="free"),
                         x = list(rot = 90)),  
           layout = c(4,1 ), 
           auto.key = list(columns = 1))
