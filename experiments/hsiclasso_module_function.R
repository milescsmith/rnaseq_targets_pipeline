hsiclasso <- function(
  expr_mat, #rows = variables; columns = observations
  metadata,
  only_hugo = FALSE,
  classifier_var,
  covars = NULL,
  print_results = FALSE,
  num_features = 20L,
  num_neighbors = 20L
){

  diffused_class <- rlang::sym(classifier_var)
  diffused_class <- rlang::enquo(diffused_class)

  expr_mat <- tibble::as_tibble(expr_mat, rownames = "symbol")

  if (!tibble::is_tibble(metadata)){
    metadata <- tibble::as_tibble(metadata, rownames = "sample_name")
  }

  if (isTRUE(only_hugo)){
    hugo_tbl <-
      tibble::as_tibble(
        x =
          HGNChelper::checkGeneSymbols(
            x = pull(expr_mat, symbol)
          )
      ) |>
      dplyr::filter(!is.na(Suggested.Symbol)) |>
      dplyr::group_by(Suggested.Symbol) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::rename(symbol = x)

    expr_mat <-
      expr_mat |>
      dplyr::filter(symbol %in% hugo_tbl[["symbol"]]) |>
      pivot_longer(
        -symbol,
        names_to = "sample_name",
        values_to = "expr"
      ) |>
      pivot_wider(
        names_from = "symbol",
        values_from = "expr"
      )
  } else {
    expr_mat <-
      expr_mat |>
      pivot_longer(
        -symbol,
        names_to = "sample_name",
        values_to = "expr"
      ) |>
      pivot_wider(
        names_from = "symbol",
        values_from = "expr"
      )
  }

  expr_data <-
    dplyr::left_join(
      x = expr_mat,
      y =
        dplyr::select(
          metadata,
          sample_name,
          {{classifier_var}}
        )
    ) |>
    dplyr::mutate(class = forcats::as_factor({{classifier_var}})) |>
    dplyr::select(-{{classifier_var}}) |>
    tibble::column_to_rownames(var = "sample_name")

  phl <- reticulate::import("pyhsiclasso", delay_load = TRUE)
  hl <- phl$HSICLasso()
  hl$input(expr_data)

  if (!is.null(covars)){
    covar_mat =
      metadata |>
      dplyr::select(tidyselect::any_of(covars)) |>
      dplyr::mutate(
        dplyr::across(
          !where(is.numeric),
          forcats::as_factor
        ),
        dplyr::across(
          where(is.factor),
          as.integer
        )
      ) |>
      as.matrix()

    hl$classification(
      num_feat      = as.integer(num_features),
      max_neighbors = as.integer(num_neighbors),
      covars        = covar_mat
    )
  } else {
    hl$classification(
      num_feat      = as.integer(num_features),
      max_neighbors = as.integer(num_neighbors)
    )
  }

  hsic_modules =
    purrr::map(
      .x = seq(num_features),
      .f = \(x)
      hl$featname[unlist(purrr::map(hl$A_neighbors[x], magrittr::add, 1))]
    ) |>
    rlang::set_names(nm = paste0("HLM", seq(num_features))) |>
    purrr::discard(
      .p = \(x) length(x) == 0
    )

  A =
    purrr::map_chr(
      .x = seq(num_features),
      .f = \(x)
      hl$featname[unlist(purrr::map(hl$A[x], magrittr::add, 1))]
    ) |>
    purrr::discard(.p = is.na)

  if (isTRUE(print_results)){
    hl$dump()
  }

  list(
    top_features = A,
    modules = hsic_modules
  )
}


#### HERE ENDS THE FUNCTION ####

tar_load(vsc_exprs)
tar_load(study_md)



hsic_modules_long <-
  tibble::enframe(hsic_modules) |>
  tidyr::unnest(cols = value) |>
  dplyr::rename(
    module = "name",
    feature = "value"
  )

vsc_hsic <-
  vsc |>
  dplyr::filter(feature %in% hsic_modules_long$feature) |>
  dplyr::left_join(hsic_modules_long)

duplicate_features <-
  vsc_hsic |>
  dplyr::group_by(feature) |>
  dplyr::summarize(number = dplyr::n()) |>
  dplyr::filter(number != 22) |>
  dplyr::pull(feature)

vsc_hsic |>
  dplyr::filter(!feature %in% duplicate_features) |>
  dplyr::group_by(responder, module) |>
  tidyHeatmap::heatmap(
    .row = feature,
    .column = sample_name,
    .value = expr,
    .scale = "row",
    palette_value = viridis::viridis(n=3)
  ) |>
  tidyHeatmap::add_tile(module)

top_hl_genes <- map_chr(hsic_modules, magrittr::extract, 1)

vsc_hsic |>
  dplyr::filter(feature %in% top_hl_genes) |>
  dplyr::select(-module) |>
  distinct() |>
  dplyr::group_by(responder) |>
  tidyHeatmap::heatmap(
    .row = feature,
    .column = sample_name,
    .value = expr,
    .scale = "row",
    palette_value = viridis::viridis(n=3)
  )
tidyHeatmap::add_tile(module)


vsc_hsic_wide <-
  vsc_hsic |>
  dplyr::group_by(sample_name, feature) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(sample_name, feature, expr) |>
  tidyr::pivot_wider(
    names_from = "feature",
    values_from = "expr"
  )

tirosh_score_modules(
  expr_obj = vsc_exprs,
  module_list = hsic_modules
) |>
  tibble::as_tibble(rownames = "sample_name") |>
  tidyr::pivot_longer(
    cols = tidyselect::starts_with("HLM"),
    names_to = "module",
    values_to = "expr"
  ) |>
  dplyr::left_join(
    dplyr::select(
      md_with_module_scores,
      -tidyselect::matches("^M[[:digit:]]"),
      -tidyselect::matches("^ldg"),
      -tidyselect::matches("^ME"),
      -mg
    )
  ) |>
  dplyr::group_by(responder) |>
  tidyHeatmap::heatmap(
    .row = module,
    .column = sample_name,
    .value = expr,
    # .scale = "row",
    palette_value = viridis::viridis(n=3)
  )

hsic_tirosh_module_scores <-
  tirosh_score_modules(
    expr_obj = vsc_exprs,
    module_list = hsic_modules
  ) |>
  tibble::as_tibble(rownames = "sample_name") |>
  tidyr::pivot_longer(
    cols = tidyselect::starts_with("HLM"),
    names_to = "module",
    values_to = "expr"
  ) |>
  dplyr::left_join(
    dplyr::select(
      md_with_module_scores,
      -tidyselect::matches("^M[[:digit:]]"),
      -tidyselect::matches("^ldg"),
      -tidyselect::matches("^ME"),
      -mg
    )
  )

ms1 <-
  hsic_tirosh_module_scores |>
  ggplot(
    aes(
      x = responder,
      y = expr,
      color = responder
    )
  ) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(facets = vars(module), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

library(moduleScoreR)

hsic_eigengene_scores <-
  scoreEigengenes(
    vsc_exprs,
    module_list = hsic_modules
  ) |>
  pivot_longer(
    cols = starts_with("HLM"),
    names_to = "module",
    values_to = "expr"
  ) |>
  left_join(
    select(
      md_with_module_scores,
      -matches("^M[[:digit:]]"),
      -matches("^ldg"),
      -matches("^ME"),
      -mg
    ),
    by = c("sample" = "sample_name")
  )

ms2 <-
  hsic_eigengene_scores |>
  ggplot(
    aes(
      x = responder,
      y = expr,
      color = responder
    )
  ) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(facets = vars(module), scales = "free_y") +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

library(ggpubr)
library(rstatix)

hsic_eigengene_stats <-
  hsic_eigengene_scores |>
  group_by(module) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()

ggboxplot(
  data = hsic_eigengene_scores,
  x = "responder",
  y = "expr",
  fill = "responder",
  facet.by = "module"
) +
  stat_pvalue_manual(
    data = hsic_eigengene_stats,
    hide.ns = TRUE
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

hsic_tirosh_stats <-
  hsic_tirosh_module_scores |>
  group_by(module) |>
  wilcox_test(formula = expr ~ responder) |>
  add_xy_position()
add_xy_positions()

ggboxplot(
  data = hsic_tirosh_module_scores,
  x = "responder",
  y = "expr",
  fill = "responder",
  add = "jitter",
  facet.by = "module",
  scales = "free_y"
) +
  stat_pvalue_manual(
    data = hsic_tirosh_stats,
    hide.ns = TRUE
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2))

future::plan(future::multisession)

hsic_GO_enrichment <-
  furrr::future_map(
    .x = hsic_modules,
    .f = \(x)
    clusterProfiler::enrichGO(
      gene = x,
      OrgDb = "org.Hs.eg.db",
      keyType = "SYMBOL",
      ont = "ALL"
    ),
    .progress = TRUE
  )|>
  keep(\(x) nrow(as.data.frame(x)) > 0) |>
  keep(\(y) all(as.data.frame(y)[["Count"]] > 1))

hsic_reactome <-
  furrr::future_map(
    .x = hsic_modules,
    .f = \(x)
    ReactomePA::enrichPathway(
      gene          =
        clusterProfiler::bitr(
          geneID = x,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = "org.Hs.eg.db"
        )[["ENTREZID"]],
      organism      = "human",
      pAdjustMethod = "fdr",
      readable      = TRUE
    ),
    .progress = TRUE
  ) |>
  keep(\(x) nrow(as.data.frame(x)) > 0) |>
  keep(\(y) all(as.data.frame(y)[["Count"]] > 1))

#### Using all within reticulate ####
library(reticulate)
use_condaenv("reticulate")
phl <- import("pyhsiclasso")
hl <- phl$HSICLasso()

vsc <-
  vsc_exprs |>
  as_tibble(rownames = "symbol") |>
  left_join(hugo_genes_tbl) |>
  filter(!is.na(Suggested.Symbol)) |>
  group_by(Suggested.Symbol) |>
  slice_sample(n = 1) |>
  select(-symbol, -Approved) |>
  column_to_rownames(var = "Suggested.Symbol")

pvsc <-
  vsc %>%
  t() %>%
  as.data.frame() %>%
  merge(
    tibble::column_to_rownames(
      md_with_module_scores[, c("sample_name", "responder")],
      var = "sample_name"
    ),
    by = 0
  ) %>%
  dplyr::rename(class = responder) %>%
  dplyr::select(-Row.names)

rm(hl)
phl <- import("pyhsiclasso")
hl <- phl$HSICLasso()
hl$input(pvsc)
hl$classification(20, max_neighbors=20L)

hsic_modules =
  purrr::map(
    .x = seq(20),
    .f = \(x)
    hl$featname[unlist(purrr::map(hl$A_neighbors[x], magrittr::add, 1))]
  ) |>
  rlang::set_names(nm = paste0("HLM", seq(20)))

top_hl_genes <- map_chr(hsic_modules, magrittr::extract, 1)

hsic_modules_long <-
  tibble::enframe(hsic_modules) |>
  tidyr::unnest(cols = value) |>
  dplyr::rename(
    module = "name",
    feature = "value"
  )

vsc_hsic <-
  pvsc |>
  as_tibble(rownames = "feature") |>
  dplyr::filter(feature %in% hsic_modules_long$feature) |>
  dplyr::left_join(hsic_modules_long) |>
  dplyr::rename(responder = "class")

vsc_hsic |>
  dplyr::filter(feature %in% top_hl_genes) |>
  dplyr::select(-module) |>
  distinct() |>
  dplyr::group_by(responder) |>
  tidyHeatmap::heatmap(
    .row = feature,
    .column = sample_name,
    .value = expr,
    .scale = "row",
    palette_value = viridis::viridis(n=3)
  )

vsc_feature_genes <-
  vsc_hsic |>
  dplyr::filter(feature %in% top_hl_genes) |>
  dplyr::select(-module) |>
  distinct()

vsc_feature_stats <-
  vsc_feature_genes |>
  dplyr::group_by(feature) |>
  rstatix::wilcox_test(formula = expr ~ responder, p.adjust.method = "fdr") |>
  rstatix::add_y_position(step.increase = 0.05)

fold_change_values <-
  vsc_feature_genes %>%
  group_by(
    feature,
    responder
  ) |>
  summarize(avg_expr = mean(expr)) |>
  pivot_wider(
    names_from = "responder",
    values_from = "avg_expr"
  ) |>
  group_by(feature) |>
  summarize(lfc = log2(responder/non_responder)) |>
  left_join(
    select(vsc_feature_stats, feature, p)
  )

ggboxplot(
  data = vsc_feature_genes,
  x = "responder",
  y = "expr",
  fill = "responder",
  add = "jitter",
  facet.by = "feature",
  scales = "free_y"
) +
  stat_pvalue_manual(
    data = vsc_feature_stats,
    hide.ns = TRUE
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_y_continuous(expand = expansion(mult = 0.2))


