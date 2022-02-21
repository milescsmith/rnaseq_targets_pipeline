plot_gex_violins <- function(
  dataset,
  goi,
  stats_data,
  x_var,
  y_var,
  palette = "ggsci::uniform_startrek"
){

  ggviolin(
    data = filter(.data = dataset, gene == goi),
    x = x_var,
    y = y_var,
    fill = x_var,
    draw_quantiles = c(0.25, 0.5, 0.75),
    trim = TRUE,
    add.params = list(
      size = 0.75,
      alpha = 0.5
      )
    ) +
    geom_quasirandom(size = 1) +
    scale_fill_paletteer_d(palette) +
    facet_wrap(
      facets = vars(gene),
      scales = "free_y"
    ) +
    stat_pvalue_manual(
      data = filter(.data = stats_data, gene == goi),
      label = "p.adj.signif",
      hide.ns = TRUE
    ) +
    labs(
      title = goi,
      caption =
        str_glue(
          "Mann - Whitney, Benjamini - Hochberg adjusted value.",
          "<br> * <i>p</i> <0.05, ",
          "** <i>p</i> <0.01, ",
          "*** <i>p</i> <0.001, ",
          "**** <i>p</i> <0.0001"
        )
    ) +
    guides(fill = guide_legend(nrow = 3)) +
    theme(axis.text.x =
            element_text(
              angle = 45,
              hjust = 1,
              vjust = 1
              ),
          plot.caption = element_markdown())
  }


extract_pathway_exprs <- function(exprs_mat, gene_list, metadata){
  sliced_expr_tbl <-
    as_tibble(
      x = exprs_mat,
      rownames = "gene"
      ) %>%
      filter(gene %in% gene_list) %>%
      pivot_longer(
        -gene,
        names_to = "sample_name",
        values_to = "expression"
        ) %>%
      left_join(
        y =
          select(
            .data = metadata,
            sample_name,
            study_group,
            grant_defined_severity,
            project_group
            )
      )

  sliced_expr_tbl
}


calc_gex_stats <- function(expr_tbl){
  sliced_expr_tbl_stats <-
    expr_tbl %>%
    group_by(gene) %>%
    wilcox_test(
      formula = expression ~ study_group
    ) %>%
    adjust_pvalue(
      method = "BH"
    ) %>%
    add_significance() %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = expr_tbl,
      group_var = gene,
      compare_value = expression,
      percent_shift_down = 0.99
    )

  sliced_expr_tbl_stats
}

grouped_complex_heatmap <- function(
  exprs_mat,
  exprs_stats,
  sample_info_tbl,
  pathway,
  gene_list,
  group_var,
  group1,
  group2,
  group1_color = "#440154FF",
  group2_color = "#FDE725FF",
  row_name_fontsize = 12
){

  if(is.character(group_var)){
    diffused_group_var <- sym(group_var)
  }

  diffused_group_var <- enquo(diffused_group_var)

  ha_vec <-
    exprs_stats %>%
    pluck(pathway) %>%
    select(
      gene,
      p.adj.signif
    )

  ha_stats <-
    gene_list %>%
    pluck(pathway) %>%
    enframe() %>%
    rename(gene = value) %>%
    filter(
      gene %in% exprs_mat[[pathway]][["gene"]],
      !gene %in% ha_vec[["gene"]]
    ) %>%
    select(-name) %>%
    mutate(p.adj.signif = "ns") %>%
    bind_rows(ha_vec) %>%
    column_to_rownames("gene")

  group1_pathway_exprs_mat <-
    exprs_mat %>%
    pluck(pathway) %>%
    select(
      gene,
      sample_name,
      expression,
    ) %>%
    pivot_wider(
      names_from = "gene",
      values_from = "expression"
    ) %>%
    column_to_rownames("sample_name") %>%
    as.matrix() %>%
    scale() %>%
    magrittr::extract(
      filter(
        .data =
          pluck(
            .x = exprs_mat,
            pathway
          ),
        {{diffused_group_var}} == group1
      ) %>%
        pull(sample_name) %>%
        unique(),
    ) %>%
    t()

  group2_pathway_exprs_mat <-
    exprs_mat %>%
    pluck(pathway) %>%
    select(
      gene,
      sample_name,
      expression,
    ) %>%
    pivot_wider(
      names_from = "gene",
      values_from = "expression"
    ) %>%
    column_to_rownames("sample_name") %>%
    as.matrix() %>%
    scale() %>%
    magrittr::extract(
      filter(
        .data =
          pluck(
            .x = exprs_mat,
            pathway
          ),
        {{diffused_group_var}} == group2
        # study_group == group2
      ) %>%
        pull(sample_name) %>%
        unique(),
    ) %>%
    t()

  group2_heatmap <-
    Heatmap(
      matrix = group2_pathway_exprs_mat,
      right_annotation = rowAnnotation(padj = anno_text(ha_stats %>% as_tibble(rownames = "gene") %>% arrange(gene) %>% deframe(), gp = gpar(fontsize = row_name_fontsize))),
      top_annotation = columnAnnotation(group_var = anno_block(gp = gpar(fill=group2_color))),
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_split =
        filter(
          .data = sample_info_tbl,
          # {{diffused_group_var}} == group2
          study_group == group2
        ) %>%
        column_to_rownames("sample_name"),
      col = circlize::colorRamp2(seq(-4,4), viridis::viridis(length(seq(-4,4)))),
      name = pathway,
      row_names_gp = gpar(fontsize = row_name_fontsize)
    )

  group1_heatmap <-
    Heatmap(
      matrix = group1_pathway_exprs_mat,
      top_annotation = columnAnnotation(group_var = anno_block(gp = gpar(fill=group1_color))),
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = FALSE,
      column_split =
        filter(
          .data = sample_info_tbl,
          {{diffused_group_var}} == group1
          # study_group == group1
        ) %>%
        column_to_rownames("sample_name"),
      col = circlize::colorRamp2(seq(-4,4), viridis::viridis(length(seq(-4,4)))),
      name = pathway,
      row_names_gp = gpar(fontsize = row_name_fontsize)
    )

  ht_list = group1_heatmap + group2_heatmap
  m1 = group1_pathway_exprs_mat
  m2 = group2_pathway_exprs_mat
  rg = range(c(group1_pathway_exprs_mat, group2_pathway_exprs_mat))
  rg[1] = rg[1] - (rg[2] - rg[1])* 0.02
  rg[2] = rg[2] + (rg[2] - rg[1])* 0.02

  anno_multiple_boxplot = function(index) {
    nr = length(index)

    pushViewport(
      viewport(
        xscale = rg,
        yscale = c(0.5, nr + 0.5)
      )
    )

    for(i in seq_along(index)) {
      grid.rect(
        y = nr-i+1,
        height = 1,
        default.units = "native"
      )

      grid.boxplot(
        m1[ index[i], ],
        pos = nr-i+1 + 0.2,
        box_width = 0.3,
        gp = gpar(fill = group1_color),
        direction = "horizontal"
      )
      grid.boxplot(
        m2[ index[i], ],
        pos = nr-i+1 - 0.2,
        box_width = 0.3,
        gp = gpar(fill = group2_color),
        direction = "horizontal"
      )
    }
    grid.xaxis()
    popViewport()
  }

  ht_list =
    rowAnnotation(
      boxplot = anno_multiple_boxplot,
      width = unit(2, "cm"),
      show_annotation_name = FALSE
    ) +
    ht_list

  lgd <-
    Legend(
      labels = c(group1, group2),
      title = "Expression",
      legend_gp = gpar(fill = c(group1_color, group2_color))
    )

  draw(
    object = ht_list,
    padding = unit(c(20, 2, 2, 2), "mm"),
    heatmap_legend_list = list(lgd)
  )
}

empty_enrichment <- function(er){
  sig_er <- dplyr::filter(er@result, p.adjust <= er@pvalueCutoff)
  if (nrow(sig_er) > 0){
    FALSE
  } else {
    TRUE
  }
}
