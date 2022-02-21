extract_viral_expression <- function(annotations, exprs, dds){
  viral_transcripts <-
    filter(
      .data = annotations,
      str_detect(
        string = transcript,
        pattern = "^ENST",
        negate = TRUE
      )
    ) %>%
    pull(gene_name) %>%
    intersect(rownames(exprs))

  detected_viral_transcripts =
    counts(dds) %>%
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
      exprs[detected_viral_transcripts,] %>%
      t() %>%
      as_tibble(rownames = "sample_name")
}
