install("milescsmith/PLIER@1.0.0")
librarian::shelf(
  c(
    "tidyverse",
    "DESeq2",
    "PLIER",
    "rsvd",
    "paletteer",
    "viridis",
    "ggpubr",
    "rstatix",
    "targets"
  )
)

grouped_add_xy_positions <- function(stats_tbl,
                                     data_tbl,
                                     group_var,
                                     compare_value,
                                     cutoff = 0.05,
                                     step_increase = 0.1){

  unique_groups <- stats_tbl %>% pull({{group_var}}) %>% unique()

  data_min_max <-
    data_tbl %>%
    select({{group_var}}, {{compare_value}}) %>%
    group_by({{group_var}}) %>%
    summarise(max = max({{compare_value}}),
              min = min({{compare_value}}),
              span = max-min,
              step = span * step_increase)

  tbl_with_positions <- map_dfr(unique_groups, function(x){
    stats_subset <- stats_tbl %>% filter({{group_var}} == x) %>% add_x_position()

    if ("p.adj" %in% names(stats_subset)){
      stats_subset <- stats_subset %>% filter(p.adj <= cutoff)
    } else {
      stats_subset <- stats_subset %>% filter(p <= cutoff)
    }

    min_max_subset <- data_min_max %>% filter({{group_var}} == x)
    if (nrow(stats_subset) > 1){
      positions <-
        seq(
          from = min_max_subset[['max']],
          by = min_max_subset[['step']],
          to = min_max_subset[['max']] + nrow(stats_subset)*min_max_subset[['step']])
      stats_subset[['y.position']] <- positions[2:length(positions)] * 0.8
    } else {
      stats_subset[["y.position"]] <- (min_max_subset[["max"]] + nrow(stats_subset)*min_max_subset[['step']]) * 0.8
    }
    stats_subset
  })

  tbl_with_positions
}

data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)

allPaths=combinePaths(canonicalPathways)

tar_load(vsd_exprs)

plierResult <- PLIER(vsd_exprs, allPaths)

targets::tar_load(annotation_info)
targets::tar_load(group_pal)

sample_disease_order_info <-
  annotation_info %>%
  as_tibble(rownames = "subject") %>%
  arrange(responder)

sample_disease_order <-
  sample_disease_order_info %>%
  pull(subject)

sample_disease_order_count <-
  sample_disease_order_info %>%
  dplyr::count(responder) %>%
  pull(n) %>%
  accumulate(sum)

pheatmap::pheatmap(
  mat = t(plierResult$B)[sample_disease_order,],
  scale="column",
  cluster_rows = FALSE,
  annotation_row=annotation_info,
  annotation_colors=group_pal,
  gaps_row = sample_disease_order_count,angle_col = 315
  )

tar_load(dds_with_scores)

plier_modules <-
  t(plierResult$B) %>%
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
      cols = 2:15,
      names_to = "module",
      values_to = "score"
    ) %>%
  dplyr::rename(response = responder) %>%
  mutate(response = fct_recode(response, responder = "1", `non-responder`= "0"))

plier_stats <-
  plier_modules %>%
  group_by(module) %>%
  wilcox_test(
    as.formula(str_glue("score ~ response")),
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl = .,
    data_tbl = plier_modules,
    group_var = module,
    compare_value = score
  )

plier_modules %>%
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
    data = plier_stats,
    label =
      ifelse(
        test = ("p.adj.signif" %in% plier_stats),
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

plotU(plierResult, auc.cutoff = 0.6, fdr.cutoff = 0.2, top = 3)

indexToPlot <-
  which(
    apply(
      X = plierResult$Uauc*(plierResult$Up<0.001),
      MARGIN = 2,
      FUN = max) > 0.75
    )

plotTopZ(
  plierRes = plierResult,
  data = vsd_exprs,
  priorMat = allPaths,
  top = 5,
  #index = indexToPlot
  )

c5_mf <- clusterProfiler::read.gmt("references/c5.go.mf.v7.4.symbols.gmt")
c5_mf_mat <-
  c5_mf %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = "term",
    values_from = "present",
    values_fill = 0
    ) %>%
  column_to_rownames("gene") %>%
  as.matrix()

plier_c5_mf <- PLIER(vsd_exprs, c5_mf_mat)
pheatmap::pheatmap(
  mat = t(plier_c5_mf$B)[sample_disease_order,],
  scale="column",
  cluster_rows = FALSE,
  annotation_row=annotation_info,
  annotation_colors=group_pal,
  gaps_row = sample_disease_order_count,angle_col = 315
)
plier_c5_mf_modules <-
  t(plier_c5_mf$B) %>%
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
    cols = 2:15,
    names_to = "module",
    values_to = "score"
  ) %>%
  dplyr::rename(response = responder) %>%
  mutate(response = fct_recode(response, responder = "1", `non-responder`= "0"))


plier_c5_mf_stats <-
  plier_c5_mf_modules %>%
  group_by(module) %>%
  wilcox_test(
    as.formula(str_glue("score ~ response")),
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl = .,
    data_tbl = plier_c5_mf_modules,
    group_var = module,
    compare_value = score
  )

plier_c5_mf_modules %>%
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
    data = plier_c5_mf_stats,
    label =
      ifelse(
        test = ("p.adj.signif" %in% plier_c5_mf_stats),
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

LV6_top_genes <-
  plier_c5_mf %>%
  chuck("Z") %>%
  as_tibble(
    rownames = "gene",
    .name_repair = \(x) as.character(paste0("LV", seq(x)))
  ) %>%
  top_n(40, LV6) %>%
  pull(gene)

LV6_data <-
  counts(dds_with_scores, normalized = TRUE) %>%
  as_tibble(rownames = "gene") %>%
  filter(gene %in% LV6_top_genes) %>%
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

LV6_stats <-
  LV6_data %>%
  group_by(gene) %>%
  wilcox_test(
    as.formula(str_glue("expr ~ response")),
    p.adjust.method = "BH"
  ) %>%
  grouped_add_xy_positions(
    stats_tbl = .,
    data_tbl = LV6_data,
    group_var = gene,
    compare_value = expr
  )

LV6_data %>%
  ggviolin(
    x = "response",
    y = "expr",
    fill = "response",
    draw_quantiles = c(0.25, 0.5, 0.75),
    add = "jitter",
    trim = TRUE,
    add.params = list(size = 1, alpha = 0.25)
  ) +
  scale_fill_manual(values = c("#CC2C00FF", "#5C88DAFF")) +
  stat_pvalue_manual(
    data = LV6_stats,
    label =
      ifelse(
        test = ("p.adj.signif" %in% LV6_stats),
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
