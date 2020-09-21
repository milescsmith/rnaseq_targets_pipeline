loadd(dds_with_scores, study_md)

cytokines <-
  colData(dds_with_scores) %>%
  as_tibble(rownames = "sample_name") %>%
  left_join(select(study_md, sample_name, cluster)) %>%
  select(
    sample_name,
    disease_class,
    cluster,
    one_of(cytokine_names),
    `BCL/CXCL13` = `BCL.CXCL13`,
    `GM-CSF` = `GM.CSF`
    )

cytokine_tbl <-
  cytokines %>%
  mutate(
    disease_class =
      fct_recode(
        .f = disease_class,
        CLE = "SCLE", CLE = "DLE", CLE = "TLE", CLE = "LP"
        )
    ) %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "cytokine"
  ) %>%
  mutate(
    value = replace_na(data = value, replace = 0),
  ) %>%
  filter(sample_name %in% study_md$sample_name)

cytokine_cluster_tbl_summary <-
  cytokine_tbl %>%
  group_by(cluster, cytokine) %>%
  summarise(
    avg = mean(value),
    iqr = IQR(value),
    n = n()
    ) %>%
  distinct()

cytokine_cluster_tbl_formatted_tbl <-
  cytokine_cluster_tbl_summary %>%
  mutate(
    desired =
      str_glue("{round(avg, 1)}\n({round(avg-iqr, 1)}-{round(avg+iqr, digits = 1)})") %>%
      as.character(),
    cluster = str_glue("Cluster {cluster}\n(N={n})")
  ) %>%
  select(
    cluster,
    cytokine,
    desired
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = "cluster",
    values_from = "desired"
  )

cytokine_cluster_pvals <-
  map_dfr(unique(cytokine_tbl$cytokine), function(i){
    cytokine_tbl %>%
      filter(cytokine == i) %>%
      kruskal_test(value ~ cluster) %>%
      mutate(
        cytokine = i,
        p = round(p, 2))
  }) %>%
  select(
    cytokine,
    `P-value` = p
    )

cytokine_cluster_compiled_tbl <-
  left_join(
    cytokine_cluster_tbl_formatted_tbl,
    cytokine_cluster_pvals
  ) %>%
  rename(`Soluble Mediator median (IQR)` = cytokine)

cytokine_disease_tbl_summary <-
  cytokine_tbl %>%
  group_by(disease_class, cytokine) %>%
  summarise(
    avg = mean(value),
    iqr = IQR(value),
    n = n()
  ) %>%
  distinct()

cytokine_disease_tbl_formatted_tbl <-
  cytokine_disease_tbl_summary %>%
  mutate(
    desired =
      str_glue("{round(avg, 1)}\n({round(avg-iqr, 1)}-{round(avg+iqr, digits = 1)})") %>%
      as.character()
  ) %>%
  select(
    disease_class,
    cytokine,
    desired
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = "disease_class",
    values_from = "desired"
  )

cytokine_disease_pvals <-
  map_dfr(unique(cytokine_tbl$cytokine), function(i){
    cytokine_tbl %>%
      filter(cytokine == i) %>%
      wilcox_test(value ~ disease_class, p.adjust.method = "BH") %>%
      mutate(
        cytokine = i,
        p = round(p, 2))
  }) %>%
  select(
    cytokine,
    `P-value` = p
  )

cytokine_disease_compiled_tbl <-
  left_join(
    relocate(cytokine_disease_tbl_formatted_tbl, Control, .after = CLE),
    cytokine_disease_pvals
  ) %>%
  rename(`Soluble Mediator median (IQR)` = cytokine)
