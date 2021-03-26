viral_transcripts =
  annot %>%
  filter(
    !str_detect(
      string = transcript,
      pattern = "^ENST"
    )
  ) %>%
  pull(gene_name) %>%
  intersect(rownames(vsd_exprs))

detected_viral_transcripts =
  counts(dds_with_scores) %>%
  t() %>%
  as_tibble(rownames = "sample_name") %>%
  select(
    sample_name,
    one_of(viral_transcripts)
  ) %>%
  pivot_longer(
    -sample_name,
    names_to = "transcript",
    values_to = "counts"
  ) %>%
  group_by(transcript) %>%
  summarise(total_counts = sum(counts)) %>%
  filter(total_counts > 0) %>%
  pull(transcript)

viral_exprs =
  vsd_exprs[detected_viral_transcripts,] %>%
  t() %>%
  as_tibble(rownames = "sample_name")

ifn_modules =
  annotated_modules %>%
  filter(type == "Interferon") %>%
  pull(module)

inflame_modules =
  annotated_modules %>%
  filter(type == "Inflammation") %>%
  pull(module)

module_scores_with_viral =
  study_md %>%
  select(
    sample_name,
    disease_class,
    # cell_type,
    one_of(names(banchereau_modules)),
    one_of(names(ldg_modules)),
    names(metasignature_module)
  ) %>%
  inner_join(wgcna_scores) %>%
  inner_join(viral_exprs) %>%
  inner_join(clusters)
