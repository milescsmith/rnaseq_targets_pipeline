plot_dispersion_estimate.DGEList <-
  function(
    object,
    CV = FALSE,
    ...
    ){
    f    <- ifelse(CV, sqrt, I)

    if (!"tagwise.dispersion" %in% names(object)){
      object <-
        edgeR::estimateDisp(
          y      = object,
          design = design_matrix,
          robust = TRUE
          )
    }

    glm_object <- glmFit(object)

    disp_data <-
      tibble::tibble(
        count_means = matrixStats::rowMeans2(edgeR::getCounts(object)),
        dispersions  = purrr::chuck(object, "tagwise.dispersion"),
        dispersions_fit = purrr::chuck(glm_object, "dispersion"),
        gene_name   = rownames(object),
        dispOutlier     = chuck(gof(glm_object), "outlier")
      ) %>%
      mutate(
        dispersions   = if_else(
          condition   = count_means > 0,
          true        = f(dispersions),
          false       = dispersions
        ),
        outlier_shape =
          if_else(
            condition = dispOutlier,
            true      = 1,
            false     = 16
          ) %>% forcats::as_factor(),
        outlier_size  =
          if_else(
            condition = dispOutlier,
            true      = 2 * 0.45,
            false     =  0.45
          ) %>% forcats::as_factor(),
        outlier_halo  =
          if_else(
            condition = dispOutlier,
            true      = "final",
            false     = "gene-est"
          ) %>% forcats::as_factor()
      )

    ymin <-
      disp_data %>%
      dplyr::filter(
        dispersions > 0 & !is.na(dispersions)
        ) %>%
      dplyr::pull(dispersions) %>%
      min() %>%
      log10() %>%
      magrittr::subtract(0.1) %>%
      floor() %>%
      magrittr::raise_to_power(10, .)

    disp_data <-
      disp_data %>%
      mutate(
        dispEst = pmax(dispersions, ymin)
      )

    disp_plot <-
      disp_data %>%
      ggplot(
        aes(
          x = count_means,
          y = py
        )
      ) +
      geom_point() +
      geom_point(
        aes(
          x = count_means,
          y = dispersions,
          size = outlier_size,
          shape = outlier_shape,
          color = outlier_halo
        )
      ) +
      scale_x_log10() +
      scale_y_log10() +
      scale_shape_manual(values = c(1, 16)) +
      scale_size_manual(values = c(1,2)) +
      scale_color_manual(values = c(
        "dodgerblue",
        "red",
        "black"), ) +
      geom_point(
        mapping =
          aes(
            x = count_means,
            y = dispersions_fit,
            color = "fitted"
          ),
        size = 1
      ) +
      labs(
        x = "mean of normalized counts",
        y = "dispersion",
        color = ""
      ) +
      guides(
        size = "none",
        shape = "none"
      ) +
      theme_pubr() +
      theme(
        legend.justification=c(1,0),
        legend.position=c(1,0)
      )

    return(disp_plot)
  }

plot_dispersion_estimate.DGEList <-
  function(
    object,
    process_method = c("edgeR", "limma"),
    ...
  ){
    process_method = match.arg(process_method)

    if (process_method == "edgeR"){

    }
  }
