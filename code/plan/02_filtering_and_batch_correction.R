corrected_counts =
  ComBat_seq(
    counts = counts(dds_import),
    batch = fct_drop(colData(dds_import)[[batch_variable]]),
    group = colData(dds_import)[[comparison_grouping_variable]]
  ) %>%
  `storage.mode<-`("integer")

dds_import_combat =
  DESeqDataSetFromMatrix(
    countData = corrected_counts,
    colData = colData(dds_import),
    design = study_design
  )

dds_filtered =
  dds_import %>%
  magrittr::extract(rowSums(counts(.)) > 1, ) %>%
  magrittr::extract(
    grep(
      pattern = "^RNA5",
      x = rownames(.),
      invert = TRUE,
      value = TRUE
      )
    )

## Sample QC filtering
# Remove samples that have a PC1 Z-score > 3. This matches what I was doing visually, but is vastly quicker.
outlier_qc =
  remove_outliers(
    dds = dds_filtered,
    pc1_zscore_cutoff = pc1_zscore_threshold,
    pc2_zscore_cutoff = pc2_zscore_threshold
    )

dds_qc = outlier_qc$dds

pca_qc = outlier_qc$pca

removed_outliers = outlier_qc$removed

dds =
  DESeq(
    dds_qc,
    parallel = TRUE,
    BPPARAM = BPPARAM
  )

sva_res = calc_sva(dds = dds, model_design = comparison_grouping_variable, n.sva = num_sva)
dds_processed = dds
sva_graph_data = sva_res$sva
