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
          ) := padj,
        "lfc standard error" = lfcSE
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

#' comparison_heatmap
#'
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
comparison_heatmap <-
  function(
    exprs,
    comparison_name,
    utbl,
    dtbl,
    md,
    comparison_variable,
    grouping_variable
  ){
    grouping_variable <-
    comparison_variable <- rlang::enquo(comparison_variable)

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
        {{comparison_variable}} %in%
          stringr::str_split(
            string = comparison_name,
            pattern = "_vs_",
            simplify = TRUE
            )
        ) |>
      dplyr::pull(subject)

    sample_disease_order_info <-
      annotation_info |>
      tibble::as_tibble(rownames = "subject") |>
      dplyr::filter(subject %in% people) |>
      dplyr::arrange({{comparison_variable}})

    plot_data <-
      exprs |>
      magrittr::extract(comparison_genes, people) |>
      tibble::as_tibble(rownames = "gene") |>
      tidyr::pivot_longer(
        cols = tidyselect::one_of(people),
        names_to = "subject",
        values_to = "expr"
        ) |>
      dplyr::group_by(gene) |>
      dplyr::mutate(
        variance = var(expr),
        expr = scale(expr)
        ) |>
      dplyr::filter(variance > 1e-10) |>
      dplyr::ungroup() |>
      dplyr::left_join(
        y = md,
        by = "subject"
      ) |>
      dplyr::arrange({{diffused_comparison}})

    comparison_group_breaks <-
      md |>
      dplyr::count({{diffused_comparison}}) |>
      dplyr::pull(n) |>
      purrr::accumulate(sum)

    plot_data |>
      tidyheatmap::tidy_heatmap(
        rows    = subject,
        columns = gene,
        values  = expr,
        scale   = "row",
        gaps_row = {{diffused_comparison}},
        annotation_row = colnames(select(plot_data, where(is.factor))),
        colors = viridisLite::viridis(n = 15, option = "D"),
        annotation_colors = RColorBrewer::brewer.pal()
      )

      pheatmap::pheatmap(
        #scale = 'column',
        main = i,
        fontsize = 9,
        fontsize_col = 9,
        cluster_rows = FALSE,
        show_rownames = FALSE,
        annotation_row = annotation_info[sample_disease_order, ],
        annotation_colors = group_pal,
        angle_col = "315",
        breaks = zscore_range,
        color = viridisLite::viridis(n = length(zscore_range), option = "D"),
        border_color = NA,
        gaps_row = comparison_group_breaks
      )
  }
