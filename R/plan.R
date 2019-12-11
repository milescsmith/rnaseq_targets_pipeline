#For ALE06
plan = drake_plan(
  non_project_controls = import_table(file=metadata_file,
                                      bucket="memory-beta",
                                      FUN=read_excel, 
                                      sheet = "main", 
                                      col_types = c("numeric", 
                                                    rep("text",10),
                                                    "numeric",
                                                    rep("text",6),
                                                    rep("numeric", 3),
                                                    "text",
                                                    "logical"),
                                      trim_ws = TRUE,
                                      na = "n/a",
                                      .name_repair = ~ make_clean_names) %>%
    mutate(sample_name = make_clean_names(sample_name, case = "all_caps")) %>%
    filter(initial_concentration_ng_ul > initial_concentration_threshold,
           disease_class == "Control") %>% 
    #Select the portions of the metadata that are useful:
    select(sample_name,
           study_group,
           disease_class,
           project,
           run_id,
           sample_alias,
           ord,
           sex,
           race_code,
           initial_concentration_ng_ul,
           final_concentration_ng_ul,
           rin) %>%
    mutate(project = factor(project),
           study_group = factor(study_group),
           disease_class = factor(disease_class),
           run_id = factor(run_id),
           initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200)),
  
  md =
    import_table(file=metadata_file,
                 bucket="memory-beta",
                 FUN=read_excel, 
                 sheet = "main", 
                 col_types = c("numeric", 
                               rep("text",10),
                               "numeric",
                               rep("text",6),
                               rep("numeric", 3),
                               "text",
                               "logical"),
                 trim_ws = TRUE,
                 na = "n/a",
                 .name_repair = ~ make_clean_names) %>%
    mutate(sample_name = make_clean_names(sample_name, case = "all_caps")) %>%
    filter(initial_concentration_ng_ul > 1.5,
           project %in% (projects_to_include %||% unique(.data$project)),
           project %nin% projects_to_exclude,
           disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
           disease_class %nin% disease_classes_to_exclude) %>% 
    #Select the portions of the metadata that are useful:
    select(sample_name,
           study_group,
           disease_class,
           project,
           run_id,
           sample_alias,
           ord,
           sex,
           race_code,
           initial_concentration_ng_ul,
           final_concentration_ng_ul,
           rin) %>%
    mutate(project = factor(project),
           study_group = factor(study_group),
           disease_class = factor(disease_class),
           run_id = factor(run_id),
           initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200)) %>%
    bind_rows(non_project_controls) %>%
    distinct(),
  
  #Inspect the metadata:
  md_cat_data = inspect_cat(md),
  md_num_data = inspect_num(md),
  
  # Import count data
  #Raw FASTQs were processed using Salmon, which created pseudocounts.  
  #To import those, we need a named list (e.g. a dictionary) of the location
  #of all files ending in "h5" which each list member named for the sample it represents.
  
  tx_sample_names = dir(path = seq_file_directory,
                        pattern = "quant.sf.gz",
                        recursive = TRUE,
                        full.name = TRUE) %>% 
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>% 
    str_split(pattern = "/") %>% 
    # When str_split splits a string, it makes everything before the matching pattern into an element of the returned list
    # even if there is nothing before the split - you just get an empty element
    # thus, the seventh element matches '012210101_S156_L002'
    map(function(x)`[[`(x,length(x)-1)) %>%
    # strip the sequencing run part ("_S156_L002") of the name
    str_split(pattern = '_') %>% 
    map_chr(`[[`,1) %>%
    make_clean_names(case = "all_caps"),
  
  tx_files = dir(path = seq_file_directory,
                 pattern = "quant.sf.gz",
                 recursive = TRUE,
                 full.name = TRUE) %>% 
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>%
    `names<-`(tx_sample_names) %>%
    `[`(!is.na(match(names(.), md$sample_name))),
  
  #Most data structures that support row or column names cannot tolerate duplicates.  Are there any duplicate samples that need to be fixed?
  # md_dupes = md %>% filter(sample_name %in% tx_sample_names) %>% get_dupes(sample_name),
  # 
  # dedupe = deduplicate_samples(md, tx_files),
  # deduplicated_md = dedupe$md,
  # deduplicated_tx_files = dedupe$sample,
  
  # Read in count files
  # final_md = deduplicated_md[which(deduplicated_md[['sample_name']] %in% names(deduplicated_tx_files)),] %>%
  #   column_to_rownames('sample_name'),
  final_md = md[which(md[['sample_name']] %in% names(tx_files)),] %>%
    column_to_rownames('sample_name'),
  samples = tx_files[rownames(final_md)],
  counts = tximport(samples,
                    type = "salmon",
                    txIn = TRUE,
                    txOut = FALSE,
                    tx2gene = annot,
                    importer = fread),
  
  # Process the data using DESeq2
  # The `DESeq` function run several subfunctions that take 
  # care of normalization, dispersion calculations, model fitting,
  # and differential expression analysis.  This can take quite some 
  # time, especially as the study design grows more complex.
  dds_import = DESeqDataSetFromTximport(txi = counts,
                                        colData = final_md,
                                        design = study_design), 
  dds_filtered = dds_import %>%
    `[`(rowSums(counts(.)) > 1, ) %>%
    `[`(grep(pattern = "^RNA5",
             x = rownames(.),
             invert = TRUE,
             value = TRUE),),
  
  ## Sample QC filtering
  # Remove samples that have a PC1 Z-score > 3. This matches what I was doing visually, but is vastly quicker.
  outlier_qc = remove_outliers(dds = dds_filtered,
                               zscore_cutoff = pc1_zscore_threshold),
  dds_qc = outlier_qc$dds,
  pca_qc = outlier_qc$pca,
  removed_outliers = outlier_qc$removed,
  
  # I would just use the DESeq() function, but running each contitutent separately makes it easier to recover from
  # failure or to observe progress
  
  dds_processed =
    DESeq(dds_qc,
          parallel = TRUE),
  vsd = vst(dds_processed),
  vsd_exprs = assay(vsd),
  
  sample_dists = vsd_exprs %>% t() %>% parallelDist(),
  #fig.width=12, fig.height=9
  sampleDistMatrix = as.matrix(sample_dists),
  
  # Minus the noise, actual differences and similarities are now apparent.
  sample_dendrogram = sample_dists %>% hclust() %>% as.dendrogram(),
  
  # Need to dissect the code, but it turns out that Seurat has the ability
  # to build an SNN graph off of a distance matrix, which we can then pass to the Leiden
  # algorithm to identify clusters.
  
  neighbors = FindNeighbors(object = sample_dists),
  
  # Test a range of resolutions and choose the one with the highest silhouette
  # coefficient
  optimal_resolution = target(optimize_cluster_resolution(snn_graph = neighbors$snn,
                                                   dist_mat = sample_dists,
                                                   from = 0,
                                                   to = 1.5,
                                                   by = 0.01) %>%
    top_n(1, coeff) %>%
    pull(res),
    resources = list(cores = parallel::detectCores()-4, gpus = 0)),
    
  clusters = leiden(as.matrix(neighbors$snn),
         n_iterations = 10,
         resolution_parameter = optimal_resolution,
         partition_type = "RBConfigurationVertexPartition"),
  
  leiden_res = tibble(sample = rownames(neighbors$snn),
                      cluster = as_factor(clusters)),
  
  annotation_info = as.data.frame(colData(dds_processed))[,c("disease_class",
                                                             "project",
                                                             "run_id",
                                                             "sex")] %>%
    as_tibble(rownames = "sample") %>%
    inner_join(leiden_res) %>%
    mutate(cluster = as_factor(cluster),
           project = as_factor(project)) %>%
    column_to_rownames(var="sample"),
  
  # Variation can also be examined in reduced dimensional space by PCA or UMAP:
  pca_results = prcomp_irlba(x = vsd_exprs) %>%
    `[[`("rotation") %>%
    as_tibble() %>%
    mutate(sample = colnames(vsd_exprs)) %>% 
    inner_join(as_tibble(colData(dds_processed), 
                         rownames="sample")) %>%
    inner_join(leiden_res),
  
  #fig.width=12, fig.height=6
  pca_plot = pca_results %>%
    ggplot(aes(x = PC1,
               y = PC2,
               color = disease_class)) +
    geom_point(alpha = 0.3) + 
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot(),
  
  #fig.width=12, fig.height=6
  pca_plot2 = pca_results %>%
    ggplot(aes(x = PC1,
               y = PC2,
               color = run_id)) +
    geom_point(alpha = 0.3) + 
    #geom_mark_hull(aes(fill=run_id)) +
    scale_color_paletteer_d(package="ggsci",
                            palette = "category20_d3") +
    scale_fill_paletteer_d(package="ggsci",
                           palette = 
                             "category20_d3") +
    labs(color = "run_id") +
    theme_cowplot(),
  
  #fig.width=12, fig.height=6
  umap_results = umap(t(vsd_exprs),
                      n_threads = detectCores(),
                      n_sgd_threads = detectCores(),
                      verbose = TRUE,
                      n_components = 3) %>% 
    as_tibble(.name_repair = "unique") %>%
    `colnames<-`(c("umap_1",
                   "umap_2",
                   "umap_3")) %>%
    mutate(sample = colnames(vsd_exprs)) %>% 
    inner_join(as_tibble(colData(dds_processed),
                         rownames="sample")) %>%
    inner_join(leiden_res),
  
  umap_plot = umap_results %>% 
    ggplot(aes(x = umap_1,
               y = umap_2,
               color = disease_class)) +
    geom_point(alpha=0.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(color = "disease_class") +
    theme_cowplot(),
  
  #fig.width=12, fig.height=6
  umap_plot2 = umap_results %>% 
    ggplot(aes(x = umap_1, 
               y = umap_2,
               color = run_id)) +
    geom_point(alpha=0.5) +
    scale_color_paletteer_d(package="ggsci",
                            palette = "category20_d3") +
    scale_fill_paletteer_d(package="ggsci",
                           palette = "category20_d3") +
    theme_cowplot(),
  
  res = results(
    object = dds_processed,
    contrast = c(comparison_grouping_variable, 
                 experimental_group,
                 control_group),
    alpha = 0.05,
    parallel = TRUE),
  
  down_table = 
    res %>%
    as_tibble(rownames = "gene") %>%
    filter(padj <= 0.05) %>%
    mutate(log2FoldChange = -log2FoldChange) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(25, log2FoldChange) %>%
    arrange(desc(log2FoldChange)),
  
  up_table = 
    res %>%
    as_tibble(rownames = "gene") %>%
    filter(padj <= 0.05) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(25, log2FoldChange) %>%
    arrange(desc(log2FoldChange)),
  
  degs = 
    res %>% 
    as_tibble(rownames = "gene") %>%
    filter(padj < 0.05),
  
  top_up = degs %>%
    filter(log2FoldChange > 0) %>%
    top_n(100, log2FoldChange) %>%
    pull(gene),
  top_down = degs %>%
    filter(log2FoldChange < 0) %>%
    top_n(100, -log2FoldChange) %>%
    pull(gene),
  
  #fig.width=10, fig.height=9
  ISGs = intersect(c("STAT1", "ADAR", "ABCE1", "RNASEL", "TYK2", "IFNAR1",
                     "IFNB1", "STAT2", "IFNAR2", "JAK1", "SAMHD1", 
                     "SOCS1", "SOCS3", "STAT1", "ISG20", "IFITM3", "IFITM1", "IRF9", "ISG15", 
                     "IFI6", "IFIT3", "USP18", "IP6K2", "PSMB8", "IFIT1", "IRF4", "IRF5", "IRF1", 
                     "IRF3", "IRF6", "IRF8", "IRF2", "IRF7", "IFITM2", "XAF1", "IFI27", "GBP2", 
                     "RSAD2", "MX2", "MX1", "IFIT2", "IFI35", "BST2", "OAS1", "OASL", "OAS2", "OAS3", "PTPN6", "PTPN11"),
                   rownames(vsd_exprs)),
  
  module_tbl = module_list %>%
    enframe() %>%
    unnest(cols = "value") %>%
    rename(module = name,
           gene = value) %>%
    inner_join(as_tibble(module_annot,
                         rownames="module")),
  
  module_ISGs = module_tbl %>%
    dplyr::group_by(gene) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::filter(type == "Interferon",
                  gene %in% rownames(vsd_exprs)) %>%
    dplyr::select(module, gene) %>%
    arrange(module) %>%
    column_to_rownames("gene"),
  
  dds_with_scores =
    scoreEigengenes(dds_processed, module_list = module_list, score_func = 'rsvd') %>%
    scoreEigengenes(object = ., module_list = ldg_modules, score_func = 'rsvd'),
  
  annotated_modules = module_annot %>% as_tibble(rownames="module") %>% filter(type != "Undetermined"),
  
  type_pal = paletteer_d(package="ggsci",
                         palette = "category20_d3",
                         n = length(unique(annotated_modules$type))) %>%
    `names<-`(unique(annotated_modules$type)),
  
  # we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
  run_groups =
    colData(dds_with_scores) %>%
    as_tibble() %>%
    pull(run_id) %>%
    unique(),
  
  run_id_color_set =
    colorRampPalette(
      brewer.pal(12, "Paired"))(length(run_groups)) %>%
    `names<-`(run_groups),
  
  chr_pal = c("Y" = "#0000FF", "X" = "#FF0000"),
  
  sex_pal = c("Male" = "#6666FF", "Female" = "#FF0000", "unk" = "#333333"),
  
  comparison_grouping_variable_colors = c("#000000",
                                          "#999999") %>% `names<-`(c(control_group, experimental_group)),
  
  cluster_pal = 
    colorRampPalette(
      brewer.pal(12, "Set1"))(length(levels(leiden_res$cluster))) %>%
    `names<-`(levels(leiden_res$cluster)),
  
  project_pal = 
    colorRampPalette(
      brewer.pal(12, "Set1"))(length(levels(annotation_info$project))) %>%
    `names<-`(levels(annotation_info$project)),
  
  group_pal =
    list(
      comparison_grouping_variable_colors,
      run_id_color_set,
      type_pal,
      chr_pal,
      sex_pal,
      cluster_pal,
      project_pal
    ) %>% `names<-`(c(comparison_grouping_variable, "run_id", "type", "chr", "sex", "cluster", "project")),
  
  module_scores =
    colData(dds_with_scores) %>%
    as_tibble(rownames="sample") %>%
    select(sample, 
           matches("^M[[:digit:]]+\\.")) %>%
    column_to_rownames("sample"),
  
  annotated_module_scores =
    colData(dds_with_scores) %>%
    as_tibble(rownames="sample") %>%
    select(sample, 
           one_of(annotated_modules$module)) %>%
    column_to_rownames("sample"),
  
  viral_transcripts =
    annot %>%
    filter(!str_detect(string = transcript,
                       pattern = "^ENST")) %>%
    pull(gene_name) %>%
    intersect(rownames(vsd_exprs)),
  
  viral_exprs =
    vsd_exprs[viral_transcripts,] %>%
    t() %>%
    as_tibble(rownames="sample"),
  
  ifn_modules =
    annotated_modules %>%
    filter(type == "Interferon") %>%
    pull(module),
  
  inflame_modules =
    annotated_modules %>%
    filter(type == "Inflammation") %>%
    pull(module),
  
  ifn_scores =
    colData(dds_with_scores)[,ifn_modules] %>%
    as.data.frame() %>%
    as_tibble(rownames="sample"),
  
  inflammation_scores =
    colData(dds_with_scores)[,inflame_modules] %>%
    as.data.frame() %>%
    as_tibble(rownames="sample"),
  
  ldg_scores =
    colData(dds_with_scores)[,names(ldg_modules)] %>%
    as.data.frame() %>%
    as_tibble(rownames="sample"),
  
  ifn_scores_with_viral =
    inner_join(viral_exprs, ifn_scores),
  
  inflammation_scores_with_viral =
    inner_join(viral_exprs, inflammation_scores),
  
  report = rmarkdown::render(
    knitr_in("report.rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE),

  vsd_exprs %>% 
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/filtered_normalized_stabilized_expression.csv")),
  colData(dds_with_scores) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/filtered_metadata.csv")),
  counts(dds_with_scores) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/filtered_transcript_counts.csv")),
  res %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/log_fold_changes.csv")),
  
  colData(dds_with_scores) %>%
    as_tibble(rownames="sample") %>%
    select(sample, 
           one_of(module_tbl$module),
           one_of(names(ldg_modules))) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/module_scores.csv")),
  
  disp_plot = plot_dispersion_estimate(dds_with_scores),
  
  save(annotation_info,
       final_md,
       group_pal,
       ifn_scores_with_viral,
       inflammation_scores_with_viral,
       ISGs,
       pca_results,
       top_down,
       top_up,
       umap_results,
       vsd_exprs,
       file = file_out("data/for_shiny_vis.RData"))
)
