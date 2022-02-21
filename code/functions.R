`%nin%` <- purrr::compose(`!`, `%in%`)

deduplicate_samples <- function(md, samples){
  if (nrow(get_dupes(md, sample_name)) > 0){
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples)) %>%
      mutate(sample_name = make_clean_names(string = sample_name,
                                            case = "all_caps"))
    deduplicated_samples = set_names(
      x = samples,
      nm = make_clean_names(
        string = names(samples),
        case = "all_caps"
      )
    )
  } else {
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples))
    deduplicated_samples = samples
  }

  return(list(md = deduplicated_md,
              samples = deduplicated_samples))
}


fix_antibody_values <- function(i) {
  recode(
    .x = i,
    Negative = "negative",
    POSITIVE = "positive",
    `no_val` = "no_val",
    Indeterminate = "indeterminate"
  )
}


convert_nuID_to_probeID <- function(
  object,
  rename_list
){
  featureNames(object) <-
    featureNames(object) %>%
    recode(!!!rename_list)

  return(object)
}

# Adapted from Seurat::AddModuleScore, which in turn took it from Tirosh (2006)
# TODO: this should probably be moved to the {moduleScoreR} package
tirosh_score_modules <- function(
  expr_obj,
  module_list,
  breaks = 25,
  num_ctrls = 100
) {

  features <- module_list
  name <- "module"
  cluster_length <- length(x = features)

  data_avg <- Matrix::rowMeans(x = expr_obj)
  data_avg <- data_avg[order(data_avg)]
  data_cut <-
    ggplot2::cut_number(
      x = data_avg + rnorm(n = length(data_avg)) / 1e30,
      n = num_ctrls,
      labels = FALSE,
      right = FALSE
    )

  names(x = data_cut) <- names(x = data_avg)
  ctrl_use <- vector(mode = "list", length = cluster_length)

  # for each module
  for (i in seq(cluster_length)) {
    # use only the module genes that are present in our dataset
    features_use <- features[[i]][which(features[[i]] %in% rownames(expr_obj))]

    # for each module gene
    for (j in seq_along(features_use)) {
      ctrl_use[[i]] <-
        c(
          ctrl_use[[i]],
          names(
            x = sample(
              x = data_cut[which(x = data_cut == data_cut[features_use[j]])],
              size = num_ctrls,
              replace = FALSE
            )
          )
        )
    }
  }

  ctrl_use <- lapply(X = ctrl_use, FUN = unique)
  ctrl_scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl_use),
    ncol = ncol(x = expr_obj)
  )

  for (i in 1:length(ctrl_use)) {
    features_use <- ctrl_use[[i]]
    ctrl_scores[i, ] <- Matrix::colMeans(x = expr_obj[features_use, ])
  }

  features_scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster_length,
    ncol = ncol(x = expr_obj)
  )

  for (i in 1:cluster_length) {
    features_use <- features[[i]][which(features[[i]] %in% rownames(expr_obj))]
    data_use <- expr_obj[features_use, , drop = FALSE]
    features_scores[i, ] <- Matrix::colMeans(x = data_use)
  }

  features_scores_use <- features_scores - ctrl_scores
  rownames(x = features_scores_use) <- names(module_list)
  features_scores_use <- as.data.frame(x = t(x = features_scores_use))
  rownames(x = features_scores_use) <- colnames(x = expr_obj)
  return(features_scores_use)
}


plot_scale_independence <- function(fitIndices){
  g <- ggplot(
    data = fitIndices,
    mapping = aes(
      x = Power,
      y = -sign(slope) * SFT.R.sq
    )
  ) +
    geom_text(
      mapping = aes(
        label = Power,
        color = "Red"
      )
    ) +
    labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale independence"
    ) +
    geom_hline(
      aes(
        yintercept = 0.8,
        color = "Red"
        )
      ) +
    theme(legend.position = "none") +
    theme_pubr()

  g
}

plot_connectivity <- function(fitIndices){
  mean.connect <-
    ggplot(
      data = sft$fitIndices,
      mapping = aes(
        x = Power,
        y = mean.k.
        )
      ) +
    geom_text(
      mapping = aes(
        label = Power,
        color = "Red"
        )
      ) +
    labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity"
      ) +
    theme(legend.position = "none") +
    theme_pubr()

    mean.connect
}


fix_na_y_pos <- function(dataset, stats_data, y_var = "transcript_id", value_var = "tx_expression"){
  if ("y.position" %in% colnames(stats_data)){
    if (any(is.na(stats_data[["y.position"]]))){
      if(is.character(y_var)){
        y_var <- sym(y_var)
      }

      if(is.character(value_var)){
        value_var <- sym(value_var)
      }

      y_var     <- enquo(y_var)
      value_var <- enquo(value_var)

      non_na_stats_data <- filter(.data = stats_data, !is.na(y.position))
      na_stats_data     <- filter(.data = stats_data, is.na(y.position))

      new_stats_data <-
        filter(
          .data = dataset,
          {{y_var}} %in% pull(na_stats_data, {{y_var}})
        ) %>%
        group_by({{y_var}}) %>%
        summarize(y.position = max({{value_var}})*1.1) %>%
        right_join(
          select(
            .data = na_stats_data,
            -y.position
          )
        ) %>%
        bind_rows(non_na_stats_data)
    } else {
      new_stats_data <- stats_data
    }
  } else {
    if(is.character(y_var)){
      y_var <- sym(y_var)
    }

    if(is.character(value_var)){
      value_var <- sym(value_var)
    }

    y_var     <- enquo(y_var)
    value_var <- enquo(value_var)

    new_stats_data <-
      filter(
        .data = dataset,
        {{y_var}} %in% pull(stats_data, {{y_var}})
      ) %>%
      group_by({{y_var}}) %>%
      summarize(y.position = max({{value_var}})*1.1) %>%
      right_join(
        stats_data
      )
  }

  new_stats_data
}
