#### module_classification ####
module_scores_with_md =
  module_scores %>%
  left_join(as_tibble(annotation_info,
                      rownames="sample_name"))

module_cluster_split =
  initial_split(
    data = module_scores_with_md %>%
      mutate(disease_class = fct_drop(disease_class)) %>%
      select(
        cluster,
        one_of(names(banchereau_modules))
      ),
    prop = 0.75,
    strata = "cluster"
  )

module_cluster_train = training(module_cluster_split)

module_cluster_test = testing(module_cluster_split)

module_cluster_rf_cv =
  train(
    form = cluster ~ .,
    method = "parRF",
    data = module_cluster_train,
    trControl =
      trainControl(
        method = "repeatedcv",
        number = 10,
        repeats = 10,
        search = "grid",
        allowParallel = TRUE
      ),
    importance = TRUE
  )

module_cluster_rf_cv_varImp =
  varImp(
    object = module_cluster_rf_cv,
    scale = FALSE,
    importance = TRUE
  )

module_disease_class_split =
  initial_split(
    data = module_scores_with_md %>%
      mutate(disease_class = fct_drop(disease_class)) %>%
      select(
        disease_class,
        one_of(names(banchereau_modules))
      ),
    prop = 0.75,
    strata = "disease_class"
  )

module_disease_class_train = training(module_disease_class_split)

module_disease_class_test = testing(module_disease_class_split)

module_disease_class_rf_cv =
  train(
    form = disease_class ~ .,
    method = "parRF",
    data = module_disease_class_train,
    trControl =
      trainControl(
        method = "repeatedcv",
        number = 10,
        repeats = 10,
        search = "grid",
        allowParallel = TRUE
      ),
    importance = TRUE
  )

module_disease_class_rf_cv_varImp =
  varImp(
    object = module_disease_class_rf_cv,
    scale = FALSE,
    importance = TRUE
  )
