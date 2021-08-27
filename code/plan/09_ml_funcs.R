rf_classifier <- function(
  dataset,
  classification_var,
  ...,
  train_proportion   = 0.75,
  engine             = c("ranger", "randomForest", "spark"),
  resampling_method  = c("cv", "bootstrap"),
  number_resamples   = 25L,
  mtry_upper_range   = dials::unknown(),
  grid_search_levels = 10L,
  print              = TRUE
){

  v <- rlang::enquo(classification_var)

  if (!rlang::quo_is_symbol(v)){
    diffused_classification <- rlang::sym(classification_var)
  } else {
    diffused_classification <- classification_var
    classification_var <- as.character(substitute(classification_var))
  }
  # diffused_classification <- rlang::enquo(diffused_classification)

  # Create data frames for the two sets:
  message("Preparing data...")
  data_split <-
    dplyr::select(
      .data = dataset,
      {{diffused_classification}},
      ...
    ) %>%
    dplyr::mutate(
      {{diffused_classification}} := forcats::fct_drop({{diffused_classification}})
    ) %>%
    rsample::initial_split(
      prop   = train_proportion,
      strata = {{diffused_classification}}
    )

  train_data <- rsample::training(data_split)
  test_data  <- rsample::testing(data_split)

  design_formula = as.formula(
    paste(
      classification_var,
      "~",
      "."
    )
  )

  resampling_method <- match.arg(resampling_method)

  training_resamples <-
    switch(
      resampling_method,
      cv        =
        rsample::vfold_cv(
          data = train_data,
          v    = number_resamples,
          strata = {{diffused_classification}}
          ),
      bootstrap =
        rsample::bootstraps(
          data  = train_data,
          times = number_resamples,
          strata = {{diffused_classification}}
        )
    )

  engine <- match.arg(engine)
  importance <-
    switch(
      engine,
      ranger       = "impurity",
      randomForest = TRUE,
      sparklyr     = NULL
    )

  message("Creating recipe...")
  data_recipe <-
    recipes::recipe(design_formula, data = train_data) %>%
    recipes::step_corr(recipes::all_predictors()) %>%
    recipes::step_center(
      recipes::all_predictors(),
      -recipes::all_outcomes()
    ) %>%
    recipes::step_scale(
      recipes::all_predictors(),
      -recipes::all_outcomes()
    ) %>%
    recipes::step_zv(recipes::all_predictors())

  engine = match.arg(engine)

  message("Setting engine...")
  rand_spec <-
    parsnip::rand_forest(
      mtry       = dials::tune(),
      trees      = dials::tune(),
    ) %>%
    parsnip::set_engine(
      engine     = engine,
      importance = importance
      ) %>%
    parsnip::set_mode(mode = "classification")

  message("Creating hyperparameter search grid...")
  rand_grid <-
    dials::grid_regular(
      dials::finalize(
        dials::mtry(range = c(1, mtry_upper_range)),
        training_resamples
      ),
      dials::trees(),
      levels = grid_search_levels
    )

  message("Creating workflow...")
  rand_wf <-
    workflows::workflow() %>%
    workflows::add_model(rand_spec) %>%
    workflows::add_recipe(data_recipe)

  message("Testing hyperparameters...")
  rand_res <-
    rand_wf %>%
    tune::tune_grid(
      resamples = training_resamples,
      grid      = rand_grid
    )

  message("Extracting most accurate model parameters...")
  final_rand_wf <-
    rand_wf %>%
    tune::finalize_workflow(
      tune::select_best(
        rand_res,
        "accuracy"
      )
    )

  message("Fitting data...")
  final_rand_fit <-
    final_rand_wf %>%
    tune::last_fit(data_split)

  list(
    model          = final_rand_fit,
    parameter_grid = rand_grid,
    workflow       = final_rand_wf,
    predictions    =
      tune::collect_predictions(
        x    = final_rand_fit,
        data = testing_data
        )
  )
}


plot_rf_classifier <-
  function(
    rf_fit,
    classification_var,
    ...
  ){
    if (is.character(substitute(classification_var))){
      diffused_classification <- rlang::sym(classification_var)
      diffused_classification <- rlang::enquo(diffused_classification)
    } else {
      diffused_classification <- enquo(classification_var)
    }

    roc_plot <-
      rf_fit %>%
      tune::collect_predictions() %>%
      yardstick::roc_curve(
        {{diffused_classification}},
        ...
        ) %>%
      ggplot2::autoplot() +
      ggpubr::theme_pubr()

    var_imp_plot <-
      rf_fit %>%
      hardhat::extract_workflow() %>%
      hardhat::extract_fit_parsnip() %>%
      vip::vip() +
      ggpubr::theme_pubr()

    cowplot::plot_grid(
      roc_plot,
      var_imp_plot
      )
  }
