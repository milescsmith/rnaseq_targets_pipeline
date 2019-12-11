library(caret)
library(yardstick)
library(rsample)

library(doParallel)
registerDoParallel(cores = parallel::detectCores() - 4)

loadd(module_scores, annotation_info)

module_scores <- module_scores %>%
  as_tibble(rownames="sample_name") %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name")) %>%
  select(disease_class, starts_with("M"))

module_split <- initial_split(data = module_scores,
                              prop = 0.8,
                              strata = "disease_class")

module_train <- training(module_split)
module_test <- testing(module_split)

module_rf <- train(disease_class ~ .,
                   method = "rf",
                   data = module_train,
                   trControl = trainControl(method = "boot",
                                            sampling = "up"))

module_glm <- caret::train(disease_class ~ .,
                           method = "glm",
                           data = module_train,
                           preProcess = "scale",
                           trControl = trainControl(method = "boot",
                                                    sampling = "up"))

module_test %>%
  mutate(`Random forest` = predict(module_rf, module_test)) %>%
  conf_mat(truth = disease_class, estimate = "Random forest")

module_test %>%
  mutate(`glm` = predict(module_glm, module_test)) %>%
  conf_mat(truth = disease_class, estimate = "glm")


#----------------------WGCNA------------------------------#
module_scores <- MEs %>%
  as_tibble(rownames="sample_name") %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name")) %>%
  select(cluster, starts_with("ME"))

module_split <- initial_split(data = module_scores[,1:10],
                              prop = 0.8,
                              strata = "cluster")

module_train <- training(module_split)
module_test <- testing(module_split)

module_rf <- train(cluster ~ .,
                   method = "parRF",
                   data = module_train,
                   trControl = trainControl(method = "repeatedcv",
                                            number = 10,
                                            repeats = 10))

module_glm <- train(cluster ~ .,
                    method = "glm",
                    data = module_train,
                    trControl = trainControl(method = "repeatedcv",
                                             number = 10,
                                             repeats = 10))

module_test %>%
  mutate(`Random forest` = predict(module_rf, module_test)) %>%
  conf_mat(truth = disease_class, estimate = "Random forest")

module_test %>%
  mutate(`glm` = predict(module_glm, module_test)) %>%
  conf_mat(truth = cluster, estimate = "glm")
