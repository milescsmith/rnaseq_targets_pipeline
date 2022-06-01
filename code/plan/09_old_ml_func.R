random_forest_model <- function(
  dataset,
  classification_var,
  ...,
  train_proportion = 0.75,
  engine           = "randomForest",
  print            = TRUE
  ){
  if (is.character(classification_var)){
    diffused_classification <- rlang::sym(classification_var)
  } else if (rlang::is_symbol(classification_var)){
    diffused_classification <- classification_var
    classification_var <- rlang::as_string(classification_var)
  }

  diffused_classification <- rlang::enquo(diffused_classification)

  design_formula = as.formula(
    paste(
      classification_var,
      "~",
      "."
    )
  )

  data_split <-
    dataset %>%
    dplyr::select(
      {{diffused_classification}},
      ...
      ) %>%
    dplyr::mutate(
      {{diffused_classification}} := forcats::fct_drop({{diffused_classification}})
    # ) %>%
    rsample::initial_split(prop = train_proportion)

  data_recipe <-
    rsample::training(data_split) %>%
    recipes::recipe(design_formula) %>%
    recipes::step_corr(recipes::all_predictors()) %>%
    recipes::step_center(
      recipes::all_predictors(),
      -recipes::all_ID(),
      -recipes::all_outcomes()
      ) %>%
    recipes::step_scale(
      recipes::all_predictors(),
      -recipes::all_ID(),
      -recipes::all_outcomes()
      ) %>%
    recipes::prep()

  data_testing <-
    data_recipe %>%
    recipes::bake(object = rsample::testing(data_split))

  data_training <-
    recipes::juice(object = data_recipe)

  data_model <-
    parsnip::rand_forest(
      trees = 100,
      mode  = "classification"
      ) %>%
    parsnip::set_engine(engine = engine)

  data_workflow <-
    workflows::workflow() %>%
    workflows::add_recipe(recipe = data_recipe) %>%
    workflows::add_model(spec = data_model)

  data_fit <-
    fit(
      data_workflow,
      data = data_training
      )

  predicted_values <-
    data_fit %>%
      predict(data_testing) %>%
      dplyr::bind_cols(
        dplyr::relocate(
          data_testing,
          classification_var
        )
      )

  if (isTRUE(print)){
    predicted_values %>%
      yardstick::metrics(
        truth = classification_var,
        estimate = .pred_class
      ) %>%
      print()
  }

  list(
    model       = data_fit,
    workflow    = data_workflow,
    predictions = predicted_values
  )

}

boost_model <- function(
  dataset,
  classification_var,
  ...,
  train_proportion = 0.75,
  engine           = "xgboost",
  print            = TRUE
){
  if (is.character(classification_var)){
    diffused_classification <- rlang::sym(classification_var)
  } else if (rlang::is_symbol(classification_var)){
    diffused_classification <- classification_var
    classification_var <- rlang::as_string(classification_var)
  }

  diffused_classification <- rlang::enquo(diffused_classification)

  design_formula = as.formula(
    paste(
      classification_var,
      "~",
      "."
    )
  )

  data_split <-
    dataset %>%
    dplyr::select(
      {{diffused_classification}},
      ...
    ) %>%
    dplyr::mutate(
      {{diffused_classification}} := forcats::fct_drop({{diffused_classification}})
    ) %>%
    rsample::initial_split(prop = train_proportion)

  data_recipe <-
    rsample::training(x = data_split) %>%
    recipes::recipe(design_formula)%>%
    recipes::step_corr(recipes::all_predictors()) %>%
    recipes::step_center(
      recipes::all_predictors(),
      -recipes::all_outcomes()
    ) %>%
    recipes::step_scale(
      recipes::all_predictors(),
      -recipes::all_outcomes()
    ) %>%
    recipes::prep()

  data_testing <-
    data_recipe %>%
    recipes::bake(rsample::testing(data_split))

  data_training <-
    recipes::juice(object = data_recipe)

  data_model <-
    parsnip::boost_tree(
      trees = 1000,
      mode  = "classification",
      importance = "impurity"
    ) %>%
    parsnip::set_engine(engine = engine)

  data_workflow <-
    workflows::workflow() %>%
    workflows::add_recipe(recipe = data_recipe) %>%
    workflows::add_model(spec = data_model)

  data_fit <-
    fit(
      data_workflow,
      data = data_training
    )

  predicted_values <-
    data_fit %>%
    predict(data_testing) %>%
    dplyr::bind_cols(
      dplyr::relocate(
        data_testing,
        classification_var
      )
    )

  if (isTRUE(print)){
    predicted_values %>%
      yardstick::metrics(
        truth = classification_var,
        estimate = .pred_class
      ) %>%
      print()
  }

  list(
    model = data_fit,
    workflow = data_workflow,
    predictions = predicted_values
  )

}
