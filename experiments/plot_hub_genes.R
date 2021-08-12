librarian::shelf(
  tidyverse,
  DESeq2,
  targets,
  paletteer,
  ggforce,
  ggpubr,
  rstatix
)

hg <-
  vsd_exprs[deframe(wgcna_hub_genes),] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(
    cols = starts_with("AA"),
    names_to = "subject_id",
    values_to = "expr"
    ) %>%
  left_join(
    as_tibble(
      colData(dds_with_scores),
      rownames = "subject_id"
      )
    ) %>%
  select(
    -starts_with("M"),
    -starts_with("ldg")
    ) %>%
  ggviolin(
    x = "responder",
    y = "expr",
    fill = "responder",
    draw_quantiles = c(0.25, 0.5, 0.75),
    add = "jitter",
    trim = TRUE,
    add.params = list(size = 1, alpha = 0.25)
  ) +
  scale_fill_paletteer_d("ggsci::uniform_startrek", direction = -1)

hg +
  facet_wrap_paginate(
    ~ gene,
    ncol = 6,
    nrow = 6,
    page = 1,
    scales = "free_y"
  )

hg +
  facet_wrap_paginate(
    ~ gene,
    ncol = 6,
    nrow = 6,
    page = 2,
    scales = "free_y"
  )
