process_deg_kable <-
  function(
    deg_table,
    deg_table_name,
    .color = "red",
    .direction = c("down", "up")){
        .direction = match.arg(.direction)

        direction =
            switch(
                EXPR = .direction,
                down = "downregulated",
                up   = "upregulated"
            )

      deg_table %>%
      dplyr::select(-baseMean) %>%
      dplyr::filter(log2FoldChange > 0) %>%
      dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
      dplyr::mutate(
        gene =
          kableExtra::cell_spec(
            x          = gene,
            bold       = TRUE
            ),
        padj =
          kableExtra::cell_spec(
            x          = padj,
            color      = "white",
            bold       = ifelse(
              test     = pvalue < 0.05,
              yes      = T,
              no       = F
            ),
            background = kableExtra::spec_color(padj)
            ),
        log2FoldChange = formattable::color_bar(.color)(log2FoldChange)
        ) %>%
      dplyr::ungroup() %>%
      dplyr::rename(
        !!paste0(
          "pvalue",
          kableExtra::footnote_marker_symbol(1)
          ) := pvalue,
        !!paste0(
          "padj",
          kableExtra::footnote_marker_symbol(2)
          ) := padj
        ) %>%
      knitr::kable(
        escape = FALSE,
        caption =
          paste0(
            "<b>",
            deg_table_name,
            ": ",
            direction,
            " in ",
            stringr::str_split(
              string   = deg_table_name,
              pattern  = "_vs_",
              simplify = TRUE
              ) %>%
              magrittr::extract(1),
            "</b>"
            )
        ) %>%
      kableExtra::kable_styling(
        bootstrap_options =
          c(
            "striped",
            "hover",
            "condensed",
            "responsive"
            ),
        full_width = TRUE,
        fixed_thead = TRUE) %>%
      kableExtra::footnote(
        symbol =
          c(
            "Wald test p-values",
            "Benjamini–Hochberg adjusted value"
            )
        ) %>%
      kableExtra::column_spec(
        column = 2,
        color = "black"
        )
}

process_deg_flextable <-
  function(
    deg_table,
    deg_table_name
    ){
    deg_table %>%
        dplyr::select(-baseMean) %>%
        dplyr::filter(log2FoldChange > 0) %>%
        dplyr::arrange(desc(log2FoldChange)) %>%
        dplyr::ungroup() %>%
        dplyr::rename(
          Gene = gene,
          `Log-fold change\nstandard error` = lfcSE,
          `P-value` = pvalue,
          `Adjusted\nP-value` = padj) %>%
        flextable::qflextable() %>%
        flextable::footnote(
          i = 1,
          j = 4:5,
          part = "header",
          value =
            flextable::as_paragraph(
              c(
                "Wald test p-values",
                "Benjamini–Hochberg adjusted value"
                )
              ),
          ref_symbols = c("*","†")
          ) %>%
        flextable::bold(part = "header") %>%
        flextable::set_formatter(
          log2FoldChange = \(x) sprintf("%.02g", x)
          ) %>%
        flextable::compose(
          i = 1,
          j = 2,
          part = "header",
          value = flextable::as_paragraph(
            "Log",
            flextable::as_sub("2"),
            "-fold change"
            )
          ) %>%
        flextable::add_header_lines(values = deg_table_name) %>%
        flextable::set_table_properties(layout = "autofit")
}

#' @title comparisonHeatMap
#'
#' @param comparison for the given comparison_variable,
#' @param exprs numeric matrix containing gene-by-subject expression values
#' @param comparison_name character vector in the form of "case_vs_control"
#' @param utbl tibble of upregulated genes. At a minimum, has a column named "gene"
#' @param dtbl tibble of downregulated genes. At a minimum, has a column named "gene"
#' @param md tibble of metadata.  Must have a column with "subject" and one matching \code{`comparison_name`}
#' @param comparison_variable Unquoted and diffused variable that is being used for comparison, i.e. "disease_state"
#'
#' @return
#' @export
#'
#' @examples
comparisonHeatMap <- function(
  comparison,
  utbl,
  dtbl,
  expr_mat,
  md,
  comparison_variable,
  grouping_variable = NULL,
  annotation_palettes,
  ...
){
  diffused_comparison <- rlang::sym(comparison_variable)
  diffused_comparison <- rlang::enquo(diffused_comparison)

  comparison_genes <-
    dplyr::pull(
      dplyr::bind_rows(
        utbl,
        dtbl
      ),
      gene
    )

  people <-
    dplyr::filter(
      .data = md,
      {{diffused_comparison}} %in% unlist(stringr::str_split(comparison, " - ", n = 2))
      ) |>
    dplyr::pull(sample_name)

  plot_data <-
    expr_mat |>
    magrittr::extract(comparison_genes, people) |>
    tibble::as_tibble(rownames = "gene") |>
    tidyr::pivot_longer(
      cols = tidyselect::one_of(people),
      names_to  = "sample_name",
      values_to = "exprs"
        ) |>
    dplyr::group_by(gene) |>
    dplyr::mutate(
      scaled_expr = as.vector(scale(exprs)),
      variance = var(exprs)
      ) |>
    dplyr::ungroup() |>
    dplyr::filter(variance > 1e-10) |>
    dplyr::left_join(y = md)

  # TODO: rewrite for ComplexHeatmap
  if (is.null(grouping_variable)){
    tidyheatmap::tidy_heatmap(
      df = plot_data,
      main =
        stringr::str_replace(
          string      = comparison,
          pattern     = "_",
          replacement = ": "
        ) |>
        stringr::str_replace_all(
          pattern     = "_",
          replacement = " "
        ),
      rows              = sample_name,
      columns           = gene,
      values            = scaled_expr,
      fontsize          = 9,
      show_rownames     = FALSE,
      cluster_cols      = TRUE,
      cluster_rows      = TRUE,
      angle_col         = 45,
      colors            = viridisLite::viridis(n = 100, option = "A"),
      annotation_row    = colnames(dplyr::select(plot_data, where(is.factor))),
      annotation_colors = annotation_palettes,
      ...
    )
  } else {
    diffused_grouping <- rlang::sym(grouping_variable)
    diffused_grouping <- rlang::enquo(diffused_grouping)

    tidyheatmap::tidy_heatmap(
      df = dplyr::arrange(.data = plot_data, {{diffused_grouping}}),
      main =
        stringr::str_replace(
          string      = comparison,
          pattern     = "_",
          replacement = ": "
        ) |>
        stringr::str_replace_all(
          pattern     = "_",
          replacement = " "
        ),
      rows              = sample_name,
      columns           = gene,
      values            = scaled_expr,
      fontsize          = 9,
      show_rownames     = FALSE,
      gaps_row          = grouping_variable,
      cluster_cols      = TRUE,
      angle_col         = 45,
      colors            = viridisLite::viridis(n = 100, option = "A"),
      annotation_row    = colnames(dplyr::select(plot_data, where(is.factor))),
      annotation_colors = annotation_palettes,
      ...
    )
    }
}

#' @title groupedHeatMap
#'
#' @param expr_mat numeric matrix containing gene-by-subject expression values
#' @param gene_list a character vector of genes to plot
#' @param md a metadata tibble to use to annotate the heatmap. Must have a
#' column with "subject"
#' @param row_grouping column in the metadata tibble to use to arrange the
#' rows of the heatmap.  If `NULL`, then hierarchical clustering is used.
#' @param col_grouping
#' @param col_annotation
#' @param annotation_palettes a list of named lists as created by
#' \code{generatePalettes}

#' @param comparison_variable Unquoted and diffused variable that is being used for comparison, i.e. "disease_state"
#'
#' @return
#' @export
#'
#' @examples
groupedHeatMap <- function(
  expr_mat,
  gene_list,
  md,
  annotation_palettes,
  row_grouping        = NULL,
  row_annotation      = NULL,
  col_grouping        = NULL,
  col_annotation      = NULL,
  scale_exprs         = TRUE,
  ...
){

  #stupid?  probably, but the if_else block below is already
  #complicated enough
  #Also, we cannot rely on passing a NULL value as apparently
  #what knitr::knit_expand passes is ""
  if (!is.null(row_grouping)){
    if (row_grouping == "hierarchical clustering" | row_grouping == "none"){
      row_grouping <- NULL
    }
  }

  if (!is.null(col_grouping)){
    if (col_grouping == "hierarchical clustering"| col_grouping == "none"){
      col_grouping <- NULL
    }
  }

  if (is.null(row_annotation)){
    row_annotation <- row_grouping
  }

  plot_data <-
    expr_mat |>
    tibble::as_tibble(rownames = "obsv") |>
    dplyr::filter(obsv %in% gene_list) |>
    tidyr::pivot_longer(
      -obsv,
      names_to = "sample_name",
      values_to = "exprs"
    ) |>
    dplyr::group_by(obsv) |>
    dplyr::mutate(
      scaled_expr =
        dplyr::case_when(
          scale_exprs ~ as.vector(scale(exprs)),
          TRUE ~ exprs
          ),
      variance = var(exprs)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(variance > 1e-10) |>
    dplyr::left_join(y = md) |>
    conditional_left_join(
      .y = col_annotation,
      .condition = tibble::is_tibble(col_annotation),
      .negate = FALSE
      )

  if(!is.null(col_annotation)){
    if (tibble::is_tibble(col_annotation)){
      col_annotation <- colnames(col_annotation)[colnames(col_annotation) != "obsv"]
    } else if (col_annotation == "none") {
     col_annotation <- NULL
    }
  }

  if (!is.null(col_grouping) & !is.null(row_grouping)){
    diffused_cols <- rlang::sym(col_grouping)
    diffused_cols <- rlang::enquo(diffused_cols)

    diffused_rows <- rlang::sym(row_grouping)
    diffused_rows <- rlang::enquo(diffused_rows)

    plot_data <- dplyr::arrange(.data = plot_data, {{diffused_cols}}, {{diffused_rows}})

  } else if (!is.null(col_grouping)){

    diffused_cols <- rlang::sym(col_grouping)
    diffused_cols <- rlang::enquo(diffused_cols)

    plot_data <- dplyr::arrange(.data = plot_data, {{diffused_cols}})

  } else if (!is.null(row_grouping)){
    diffused_rows <- rlang::sym(row_grouping)
    diffused_rows <- rlang::enquo(diffused_rows)

    plot_data <- dplyr::arrange(.data = plot_data, {{diffused_rows}})
  }

  tidyheatmap::tidy_heatmap(
    df                = plot_data,
    rows              = sample_name,
    columns           = obsv,
    values            = exprs,
    fontsize          = 9,
    show_rownames     = FALSE,
    cluster_cols      = ifelse(test = is.null(col_grouping), yes = TRUE, no = FALSE),
    gaps_col          = col_grouping,
    cluster_rows      = ifelse(test = is.null(row_grouping), yes = TRUE, no = FALSE),
    gaps_row          = row_grouping,
    angle_col         = 45,
    colors            = viridisLite::viridis(n = 100, option = "A"),
    annotation_row    = row_annotation,
    annotation_col    = col_annotation,
    annotation_colors = annotation_palettes,
    ...
  )

}

DimPlotHull <- function(
  .data,
  grouping,
  reduction = c("pca", "umap"),
  dims = c(1,2)
){
  diffused_grouping = rlang::sym(grouping)
  diffused_grouping = rlang::enquo(diffused_grouping)

  reduction = match.arg(reduction)
  if (reduction == "pca"){
    dim1 = rlang::sym(glue::glue("PC{dims[[1]]}"))
    dim1 = rlang::enquo(dim1)

    dim2 = rlang::sym(glue::glue("PC{dims[[2]]}"))
    dim2 = rlang::enquo(dim2)
  } else if (reduction == "umap"){
    dim1 = rlang::sym(glue::glue("umap_{dims[[1]]}"))
    dim1 = rlang::enquo(dim1)

    dim2 = rlang::sym(glue::glue("umap_{dims[[2]]}"))
    dim2 = rlang::enquo(dim2)
  }

  p1 <-
    ggplot2::ggplot(
      data = .data,
      ggplot2::aes(
        x = {{dim1}},
        y = {{dim2}},
        fill = {{diffused_grouping}},
        line = 1
      ),
      shape = 1
    ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "black",
      size = 2.5,
      stroke = 0.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggpubr::theme_pubr() +
    ggplot2::theme(legend.position = "none")

  p2 <-
    ggplot2::ggplot(
      data = .data,
      ggplot2::aes(
        x = {{dim1}},
        y = {{dim2}},
        fill = {{diffused_grouping}},
        line = 1
      ), shape = 1
    ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "black",
      size = 2.5,
      stroke = 0.5,
      alpha = 0.6
    ) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggpubr::theme_pubr() +
    ggforce::geom_mark_ellipse(
      ggplot2::aes(
        fill = {{diffused_grouping}},
        color = {{diffused_grouping}}
      ),
      alpha = 0.05
    ) +
    ggplot2::theme(legend.position = "bottom")

  legend <- cowplot::get_legend(p2)

  cowplot::plot_grid(
    cowplot::plot_grid(
      p1 + ggplot2::theme(legend.position = "none"),
      p2 + ggplot2::theme(legend.position = "none")
    ),
    legend,
    nrow = 2,
    rel_heights = c(1, .2)
  )
}

violinPanel <- function(
  expr_tbl,
  grouping_col,
  stats_tbl,
  values,
  facet_var = NULL,
  var_list = NULL,
  filter_var = NULL
){

  facet_var <- rlang::sym(facet_var)
  facet_var <- rlang::enquo(facet_var)

  stats_tbl <- pluck(stats_tbl, grouping_col)

  if (!is.null(filter_var)){
    filter_var <- rlang::sym(filter_var)
    filter_var <- rlang::enquo(filter_var)
    expr_tbl <-
      dplyr::filter(
        .data = expr_tbl,
        {{filter_var}} %in% var_list
      )
    stats_tbl <-
      dplyr::filter(
        .data = stats_tbl,
        {{filter_var}} %in% var_list
      )
  }

  ggpubr::ggviolin(
    data = expr_tbl,
    x = grouping_col,
    y = values,
    fill = grouping_col,
    draw_quantiles = c(0.25, 0.5, 0.75),
    add = "jitter",
    trim = TRUE,
    add.params = list(size = 1, alpha = 0.25)
  ) +
    paletteer::scale_fill_paletteer_d("ggsci::uniform_startrek") +
    ggpubr::stat_pvalue_manual(
      data = stats_tbl,
      label =
        ifelse(
          test = ("p.adj.signif" %in% stats_tbl),
          yes = "p.adj.signif",
          no = "p"
        ),
      hide.ns = TRUE,
      bracket.size = 0.1,
      tip.length = 0.01,
      size = 4
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 9
        ),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "none"
    ) +
    ggplot2::facet_wrap(
      facets = vars({{facet_var}}),
      scales = "free_y",
      nrow = 3,
      ncol = 3,
    )
}

groupedModuleGSEAPlots <- function(
  dataset,
  module_list
){
  module_data <-
    dplyr::filter(
      .data = dataset,
      module %in% module_list
    )

  ggplot2::ggplot(
    data = module_data,
    mapping =
      ggplot2::aes(
        y = GeneRatio,
        x = order,
        color = p.adjust
      )
  ) +
  ggplot2::geom_point(
    stat = "identity",
    size = 2
  ) +
  ggplot2::scale_x_continuous(
    breaks = module_data[["order"]],
    labels =
      stringr::str_wrap(
        paste(
          module_data[["ID"]],
          module_data[["Description"]],
          sep = ":"
          ),
        width = 40
      ),
    expand = c(0.1,0.1)
  ) +
  ggplot2::scale_color_gradient(
    low = "blue",
    high = "red"
  ) +
  ggplot2::labs(color = "adjusted\np value") +
  ggplot2::facet_wrap(
    facets = ggplot2::vars(module),
    ncol = 2,
    nrow = 3,
    scales = "free",
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 9),
    axis.text.x = ggplot2::element_text(size = 9),
    panel.background =
      ggplot2::element_rect(
        fill = "white",
        linetype = "solid"
      ),
    panel.border =
      ggplot2::element_rect(
        colour = "black",
        fill=NA,
        size=1
      ),
    panel.grid = ggplot2::element_line(color = "grey90"),
    axis.title.y = ggplot2::element_blank(),
    legend.position = "below"
    ) +
  ggplot2::coord_flip()
}

groupedComplexHeatmap <- function(
  expr_mat,
  gene_list,
  md,
  annotation_palettes,
  row_grouping        = NULL,
  row_annotation      = NULL,
  col_grouping        = NULL,
  col_annotation      = NULL,
  scale_exprs         = TRUE,
  silent              = FALSE,
  ...
){

  if (!is.null(row_grouping)){
      if (row_grouping == "hierarchical clustering" | row_grouping == "none"){
        row_grouping <- NULL
      }
    }

    if (!is.null(col_grouping)){
      if (col_grouping == "hierarchical clustering"| col_grouping == "none"){
        col_grouping <- NULL
      }
    }

    if (is.null(row_annotation)){
      row_annotation <- row_grouping
    }

    plot_data <-
      expr_mat |>
      tibble::as_tibble(rownames = "obsv") |>
      dplyr::filter(obsv %in% gene_list) |>
      tidyr::pivot_longer(
        -obsv,
        names_to = "sample_name",
        values_to = "exprs"
      ) |>
      dplyr::group_by(obsv) |>
      dplyr::mutate(
        scaled_expr =
          dplyr::case_when(
            scale_exprs ~ as.vector((exprs - mean(exprs))/sd(exprs)),
            TRUE ~ exprs
          ),
        variance = var(exprs)
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(variance > 1e-10) |>
      dplyr::left_join(y = md) |>
      conditional_left_join(
        .y = col_annotation,
        .condition = tibble::is_tibble(col_annotation),
        .negate = FALSE
      )

    if(!is.null(col_annotation)){
      if (tibble::is_tibble(col_annotation)){
        col_annotation <- colnames(col_annotation)[colnames(col_annotation) != "obsv"]
      } else if (col_annotation == "none") {
        col_annotation <- NULL
      }
    }

  row_annotation_data <-
    if (!is.null(row_annotation)){
      rad <-
        dplyr::select(
          .data = plot_data,
          sample_name,
          tidyselect::one_of(row_annotation)
        ) |>
        dplyr::distinct() |>
        tibble::column_to_rownames("sample_name")
      ComplexHeatmap::rowAnnotation(
        df = rad,
        col = annotation_palettes
        )
    } else {
      NULL
    }

  col_annotation_data <-
    if (!is.null(col_annotation)){
      cad <-
        dplyr::select(
          .data = plot_data,
          obsv,
          tidyselect::one_of(col_annotation)
        ) |>
        dplyr::distinct() |>
        tibble::column_to_rownames("obsv")
      ComplexHeatmap::columnAnnotation(
        df = cad,
        col = annotation_palettes
        )
    } else {
      NULL
    }

  plot_mat <-
    dplyr::select(
      .data = plot_data,
      obsv,
      sample_name,
      scaled_expr
      ) |>
    tidyr::pivot_wider(
      names_from = "obsv",
      values_from = "scaled_expr"
    ) |>
    tibble::column_to_rownames("sample_name")

  if (!is.null(col_grouping) & !is.null(row_grouping)){
    dro <- rlang::sym(row_grouping)
    dro <- rlang::enquo(dro)

    row_order_data <-
      dplyr::select(
        .data = plot_data,
        sample_name,
        all_of(row_grouping)
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(row_grouping)

    row_count <-
      dplyr::pull(
        row_order_data,
        {{dro}}
      )
    row_order <-
      dplyr::pull(
        .data = plot_data,
        sample_name
        ) |>
      unique()

    dco <- rlang::sym(col_grouping)
    dco <- rlang::enquo(dco)

    col_order_data <-
      dplyr::select(
        .data = plot_data,
        obsv,
        all_of(col_grouping)
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(col_grouping)

    col_count <-
      dplyr::pull(
        col_order_data,
        {{dco}}
      )
    col_order <-
      dplyr::pull(
        .data = plot_data,
        obsv
        ) |>
      unique()

    plot_mat <-
      magrittr::extract(
        plot_mat,
        row_order,
        col_order
        )

  } else if (!is.null(col_grouping)){
    dco <- rlang::sym(col_grouping)
    dco <- rlang::enquo(dco)

    col_order_data <-
      dplyr::select(
        .data = plot_data,
        obsv,
        all_of(col_grouping)
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(col_grouping)

    col_count <-
      dplyr::pull(
        col_order_data,
        {{dco}}
        )

    col_order <-
      dplyr::pull(
        .data = plot_data,
        obsv
        ) |>
      unique()

    plot_mat <-
      magrittr::extract(
        plot_mat,
        ,
        col_order
        )

    row_count = NULL
  } else if (!is.null(row_grouping)){
    dro <- rlang::sym(row_grouping)
    dro <- rlang::enquo(dro)

    row_order_data <-
      dplyr::select(
        .data = plot_data,
        sample_name,
        all_of(row_grouping)
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(row_grouping)

    row_count <-
      dplyr::pull(
        row_order_data,
        {{dro}}
        )
    row_order <-
      dplyr::pull(
        .data = row_order_data,
        sample_name
        ) |>
      unique()

    col_count = NULL

    plot_mat <- magrittr::extract(plot_mat, row_order,)
  } else {
    row_count = NULL
    col_count = NULL
  }

  color_range =
    length(
      x =
        seq(
          from = floor(min(plot_data$scaled_expr)),
          to   = ceiling(max(plot_data$scaled_expr)),
          by = 1
          )
      )

  ComplexHeatmap::Heatmap(
    matrix = as.matrix(plot_mat),
    col = viridisLite::cividis(n = color_range),
    border = TRUE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_centered = TRUE,
    row_split = row_count,
    column_split = col_count,
    left_annotation = row_annotation_data,
    top_annotation = col_annotation_data,
    name = "Expression\nZ-score",
    row_title_gp = grid::gpar(fontsize = 9),
    row_names_gp = grid::gpar(fontsize = 9),
    column_title_gp = grid::gpar(fontsize = 9),
    column_title_rot = 90,
    column_names_gp = grid::gpar(fontsize = 9),
    ...
  )
}

comparisonComplexHeatMap <- function(
  comparison,
  utbl,
  dtbl,
  expr_mat,
  md,
  comparison_variable,
  grouping_variable = NULL,
  annotation_palettes,
  row_annotation = NULL,
  col_annotation = NULL,
  scale_exprs = TRUE,
  ...
){
  diffused_comparison <- rlang::sym(comparison_variable)
  diffused_comparison <- rlang::enquo(diffused_comparison)

  comparison_genes <-
    dplyr::pull(
      dplyr::bind_rows(
        utbl,
        dtbl
      ),
      gene
    )

  samples <-
    dplyr::filter(
      .data = md,
      {{diffused_comparison}} %in% unlist(stringr::str_split(comparison, " - ", n = 2))
    ) |>
    dplyr::pull(sample_name)

  plot_data <-
    expr_mat |>
    magrittr::extract(comparison_genes, samples) |>
    tibble::as_tibble(rownames = "obsv") |>
    tidyr::pivot_longer(
      cols = tidyselect::one_of(samples),
      names_to  = "sample_name",
      values_to = "exprs"
    ) |>
    dplyr::group_by(obsv) |>
    dplyr::mutate(
      scaled_expr =
        dplyr::case_when(
          scale_exprs ~ as.vector((exprs - mean(exprs))/sd(exprs)),
          TRUE ~ exprs
        ),
      variance = var(exprs)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(variance > 1e-10) |>
    dplyr::left_join(y = md)

  row_annotation_data <-
    if (!is.null(row_annotation)){
      rad <-
        dplyr::select(
          .data = plot_data,
          sample_name,
          tidyselect::one_of(row_annotation)
        ) |>
        dplyr::distinct() |>
        tibble::column_to_rownames("sample_name")
      ComplexHeatmap::rowAnnotation(
        df = rad,
        col = annotation_palettes
      )
    } else {
      NULL
    }

  col_annotation_data <-
    if (!is.null(col_annotation)){
      cad <-
        dplyr::select(
          .data = plot_data,
          obsv,
          tidyselect::one_of(col_annotation)
        ) |>
        dplyr::distinct() |>
        column_to_rownames("obsv")
      ComplexHeatmap::columnAnnotation(
        df = cad,
        col = annotation_palettes
      )
    } else {
      NULL
    }

  plot_mat <-
    dplyr::select(
      .data = plot_data,
      obsv,
      sample_name,
      scaled_expr
    ) |>
    tidyr::pivot_wider(
      names_from = "obsv",
      values_from = "scaled_expr"
    ) |>
    tibble::column_to_rownames("sample_name") |>
    as.matrix()

  plot_title <-
    stringr::str_replace(
      string      = comparison,
      pattern     = "-",
      replacement = "vs"
    )

  color_range =
    length(
      x =
        seq(
          from = floor(min(plot_data$scaled_expr)),
          to   = ceiling(max(plot_data$scaled_expr)),
          by = 1
        )
    )

  if (grouping_variable == "hierarchical_clustering"){
    grouping_variable <- NULL
  }

  # TODO: rewrite for ComplexHeatmap
  if (!is.null(grouping_variable)){
    diffused_grouping <- rlang::sym(grouping_variable)
    diffused_grouping <- rlang::enquo(diffused_grouping)

    row_order_data <-
      dplyr::select(
        .data = plot_data,
        sample_name,
        all_of(grouping_variable)
      ) |>
      dplyr::distinct() |>
      dplyr::arrange({{diffused_grouping}})

    row_count <-
      dplyr::pull(
        row_order_data,
        {{diffused_grouping}}
      )
    row_order <-
      dplyr::pull(
        .data = row_order_data,
        sample_name
      )
    plot_mat <- magrittr::extract(plot_mat, row_order,)
  } else {
    row_count <- NULL
  }

  ComplexHeatmap::Heatmap(
    matrix                = plot_mat,
    border                = TRUE,
    col                   = viridisLite::viridis(n = color_range),
    show_row_names        = FALSE,
    column_names_rot      = 45,
    column_names_centered = TRUE,
    row_split             = row_count,
    right_annotation      = row_annotation_data,
    top_annotation     = col_annotation_data,
    name                  = "Expression\nZ-score",
    column_title          = plot_title,
    row_title_gp = grid::gpar(fontsize = 9),
    row_names_gp = grid::gpar(fontsize = 9),
    column_title_gp = grid::gpar(fontsize = 9),
    column_title_rot = 90,
    column_names_gp = grid::gpar(fontsize = 9),
    ...
  )
}
