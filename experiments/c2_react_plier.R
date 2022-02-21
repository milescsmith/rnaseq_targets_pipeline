librarian::shelf(tidyverse, clusterProfiler, pheatmap, PLIER, targets, furrr, future)
plan(multisession, workers=6)

tar_load(vsd_exprs)

c2_react <- clusterProfiler::read.gmt("references/c2.cp.reactome.v7.4.symbols.gmt")
c2_react_mat <-
  c2_react %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = "term",
    values_from = "present",
    values_fill = 0
  ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

plier_c2_react <- PLIER(vsc_exprs, c2_react_mat)
pheatmap::pheatmap(
  mat = t(plier_c2_react$B)[sample_disease_order,],
  scale="column",
  cluster_rows = FALSE,
  annotation_row=annotation_info,
  annotation_colors=group_pal,
  gaps_row = sample_disease_order_count,angle_col = 315
)
plier_c2_react_modules <-
  t(plier_c2_react$B) %>%
  as_tibble(rownames = "sample_name") %>%
  inner_join(
    as_tibble(
      colData(dds_with_scores),
      rownames = "sample_name"
    ) %>%
      select(
        disease_class,
        responder,
        sex,
        sample_name
      )
  ) %>%
  pivot_longer(
    cols = 2:13,
    names_to = "module",
    values_to = "score"
  ) %>%
  dplyr::rename(response = responder) %>%
  mutate(response = fct_recode(response, responder = "1", `non-responder`= "0"))


plier_c2_react_stats <-
  plier_c2_react_modules %>%
  group_by(module) %>%
  wilcox_test(
    as.formula(str_glue("score ~ response")),
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl = .,
    data_tbl = plier_c2_react_modules,
    group_var = module,
    compare_value = score
  )

plier_c2_react_modules %>%
  ggviolin(
    x = "response",
    y = "score",
    fill = "response",
    draw_quantiles = c(0.25, 0.5, 0.75),
    add = "jitter",
    trim = TRUE,
    add.params = list(size = 1, alpha = 0.25)
  ) +
  scale_fill_paletteer_d("ggsci::uniform_startrek") +
  stat_pvalue_manual(
    data = plier_c2_react_stats,
    label =
      ifelse(
        test = ("p.adj.signif" %in% plier_c2_react_stats),
        yes = "p.adj.signif",
        no = "p"
      ),
    hide.ns = TRUE,
    bracket.size = 0.1,
    tip.length = 0.01,
    size = 4
  ) +
  theme_pubr() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
    axis.text.y = element_text(size = 9),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(
    facets = vars(module),
    scales = "free_y",
    ncol = 4
  )

LV11_top_genes <-
  plier_c2_react %>%
  chuck("Z") %>%
  as_tibble(
    rownames = "gene",
    .name_repair = \(x) as.character(paste0("LV", seq(x)))
    ) %>%
  top_n(40, LV11) %>%
  pull(gene)

LV11_data <-
  counts(dds_with_scores, normalized = TRUE) %>%
  as_tibble(rownames = "gene") %>%
  filter(gene %in% LV11_top_genes) %>%
  pivot_longer(
    cols = starts_with("AA"),
    names_to = "sample_name",
    values_to = "expr"
  ) %>%
  left_join(
    as_tibble(
      annotation_info,
      rownames="sample_name"
    )
  ) %>%
  dplyr::rename(response = responder) %>%
  mutate(response = fct_recode(response, responder = "1", `non-responder`= "0"))

LV11_stats <-
  LV11_data %>%
  group_by(gene) %>%
  wilcox_test(
    as.formula(str_glue("expr ~ response")),
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl = .,
    data_tbl = LV11_data,
    group_var = gene,
    compare_value = expr
  )

LV11_data %>%
  ggviolin(
    x = "response",
    y = "expr",
    fill = "response",
    draw_quantiles = c(0.25, 0.5, 0.75),
    add = "jitter",
    trim = TRUE,
    add.params = list(size = 1, alpha = 0.25)
  ) +
  scale_fill_paletteer_d("ggsci::uniform_startrek") +
  stat_pvalue_manual(
    data = LV11_stats,
    label =
      ifelse(
        test = ("p.adj.signif" %in% LV11_stats),
        yes = "p.adj.signif",
        no = "p"
      ),
    hide.ns = TRUE,
    bracket.size = 0.1,
    tip.length = 0.01,
    size = 4
  ) +
  theme_pubr() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
    axis.text.y = element_text(size = 9),
    strip.text.x = element_text(size = 9)
  ) +
  facet_wrap(~gene, scales = "free_y")
