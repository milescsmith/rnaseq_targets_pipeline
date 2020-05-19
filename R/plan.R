#For ALE06
analysis_plan = drake_plan(

  metadata = read_excel(path = metadata_file,
                        sheet = "main",
                        col_types = c("numeric",
                                      rep("text",10),
                                      "numeric",
                                      rep("text",6),
                                      rep("numeric", 3),
                                      "text",
                                      "text"),
                        trim_ws = TRUE,
                        na = "n/a") %>%
    clean_names() %>%
    mutate(sample_name = make_clean_names(sample_name, case = "all_caps"),
           disease_class = str_remove(
             string = recode(disease_class, "Unaffected Control" = "Control"),
             pattern = " Patient")
           ) %>%
    filter(initial_concentration_ng_ul > initial_concentration_threshold),

  non_project_controls =
    metadata %>%
    filter(disease_class == "Control") %>%
    #Select the portions of the metadata that are useful:
    select(sample_name,
           disease_class,
           project,
           run_id,
           sample_alias,
           sex,
           race_code,
           initial_concentration_ng_ul,
           final_concentration_ng_ul,
           rin),

  md =
    metadata %>%
    filter(project %in% (projects_to_include %||% unique(.data$project)),
           project %nin% projects_to_exclude,
           disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
           disease_class %nin% disease_classes_to_exclude) %>%
    #Select the portions of the metadata that are useful:
    select(sample_name,
           disease_class,
           project,
           run_id,
           sample_alias,
           sex,
           race_code,
           initial_concentration_ng_ul,
           final_concentration_ng_ul,
           rin) %>%
    bind_rows(non_project_controls) %>%
    mutate(project = as_factor(project),
           disease_class = as_factor(disease_class),
           run_id = factor(run_id),
           initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200),
           disease_class = as_factor(disease_class),
           sex = as_factor(sex)) %>%
    distinct(),

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
    set_names(tx_sample_names) %>%
    `[`(!is.na(match(names(.), md$sample_name))),

  final_md = md[which(md[['sample_name']] %in% names(tx_files)),] %>%
    mutate(disease_class = fct_relevel(disease_class, {{control_group}})) %>%
    column_to_rownames('sample_name'),

  #Inspect the metadata:
  md_cat_data = inspectdf::inspect_cat(final_md),
  md_num_data = inspectdf::inspect_num(final_md),

  samples = tx_files[rownames(final_md)],
  counts = tximport(samples,
                    type = "salmon",
                    txIn = TRUE,
                    txOut = FALSE,
                    tx2gene = annot,
                    importer = data.table::fread),

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
                               pc1_zscore_cutoff = pc1_zscore_threshold,
                               pc2_zscore_cutoff = pc2_zscore_threshold),
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

  sample_dists = vsd_exprs %>% t() %>% parallelDist::parallelDist(),
  #fig.width=12, fig.height=9
  sampleDistMatrix = as.matrix(sample_dists),

  # Minus the noise, actual differences and similarities are now apparent.
  sample_dendrogram = sample_dists %>% hclust() %>% as.dendrogram(),

  # Need to dissect the code, but it turns out that Seurat has the ability
  # to build an SNN graph off of a distance matrix, which we can then pass to the Leiden
  # algorithm to identify clusters.

  neighbors = Seurat::FindNeighbors(object = sample_dists, nn.method = "annoy"),

  # Test a range of resolutions and choose the one with the highest silhouette
  # coefficient
  cluster_res_scan = optimize_cluster_resolution(
    snn_graph = as.matrix(neighbors$snn),
    dist_mat = sample_dists,
    from = 0.1,
    to = 1.5,
    by = 0.1),

  optimal_resolution =
    cluster_res_scan %>%
    top_n(1, coeff) %>%
    pull(res),

   clusters = leiden(as.matrix(neighbors$snn),
          n_iterations = 10,
          resolution_parameter = sample(optimal_resolution, 1), # On the chance that we have multiple equally good resolutions
          partition_type = "RBConfigurationVertexPartition") %>%
    set_names(rownames(neighbors$snn)) %>%
    enframe(name = "sample_name", value = "cluster") %>%
    mutate(cluster = as_factor(cluster)),

  #clusters = fastkmed(sampleDistMatrix,
  #                    ncluster = optimal_resolution)[['cluster']] %>%
  #  enframe(name = "sample_name",
  #          value = "cluster") %>%
  #  mutate(cluster = as_factor(cluster)),

  # leiden_res = tibble(sample_name = rownames(neighbors$snn),
  #                     cluster = as_factor(clusters)),

  annotation_info = as.data.frame(colData(dds_with_scores))[,c("disease_class",
                                                             "sex")] %>%
    as_tibble(rownames = "sample_name") %>%
    left_join(clusters) %>%
    column_to_rownames(var="sample_name"),

  # Variation can also be examined in reduced dimensional space by PCA or UMAP:
  pca_results = irlba::prcomp_irlba(x = vsd_exprs) %>%
    `[[`("rotation") %>%
    as_tibble() %>%
    mutate(sample_name = colnames(vsd_exprs)) %>%
    inner_join(as_tibble(colData(dds_processed),
                         rownames="sample_name")) %>%
    inner_join(clusters),

  umap_results = uwot::umap(t(vsd_exprs),
                      n_threads = detectCores(),
                      n_sgd_threads = detectCores(),
                      verbose = TRUE,
                      n_components = 3) %>%
    as_tibble(.name_repair = "unique") %>%
    set_names(c("umap_1",
                   "umap_2",
                   "umap_3")) %>%
    mutate(sample_name = colnames(vsd_exprs)) %>%
    inner_join(as_tibble(colData(dds_processed),
                         rownames="sample_name")) %>%
    inner_join(clusters),

  disease_classes = final_md %>% filter(disease_class != "Control") %>% pull(disease_class) %>% as.character() %>% unique(),

  # res = map(disease_classes, function(i){
  #   results(
  #     object = dds_processed,
  #     contrast = c("disease_class",
  #                  i,
  #                  "Control"),
  #     alpha = 0.05,
  #     parallel = TRUE)
  #   }) %>%
  #   set_names(map_chr(disease_classes, function(i){ paste0(i, "_vs_Control")})),

  res = map(seq(from = 2, to = length(resultsNames(dds_processed))), function(i){
    lfcShrink(dds_processed,
              coef = resultsNames(dds_processed)[[i]],
              parallel = TRUE,
              type = "apeglm")
  }) %>%
  set_names(map_chr(disease_classes, function(i){ paste0(i, "_vs_Control")})),

  down_tables = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(rownames = "gene") %>%
      filter(!is.na(padj) & padj <= 0.05) %>%
      mutate(log2FoldChange = -log2FoldChange) %>%
      mutate_at(vars(-gene), list(~signif(., 2))) %>%
      top_n(25, log2FoldChange) %>%
      arrange(desc(log2FoldChange))
  }) %>%
    set_names(map_chr(disease_classes, function(i){ paste0(i, "_vs_Control")})),

  up_tables = map(seq_along(res), function(i){
    res[[i]] %>%
    as_tibble(rownames = "gene") %>%
    filter(!is.na(padj) & padj <= 0.05) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(25, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
  }) %>%
    set_names(map_chr(disease_classes, function(i){ paste0(i, "_vs_Control")})),

  degs = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(rownames = "gene") %>%
      filter(padj < 0.05) %>%
      filter(abs(log2FoldChange) >= 0.5) %>%
      pull(gene)
  }) %>%
    set_names(map_chr(disease_classes,
                      function(i){ paste0(i, "_vs_Control") })),

  deg_class = enframe(degs) %>%
    unnest(cols = c(value)) %>%
    group_by(value) %>%
    mutate(count = n()) %>%
    mutate(name = case_when(count == 1 ~ name,
                            count > 1 ~ "multiple")) %>%
    select(-count) %>%
    distinct() %>%
    column_to_rownames(var = "value") %>%
    set_names("comparison"),

  deg_means = vsd_exprs %>%
    t() %>%
    as_tibble(rownames = "name") %>%
    select(name, one_of(rownames(deg_class))) %>%
    pivot_longer(-name,
                 names_to = "gene",
                 values_to = "expr") %>%
    left_join(as_tibble(colData(dds_processed), rownames = "name") %>%
                select(name, sex, race_code, disease_class)) %>%
    group_by(gene,
             disease_class) %>%
    summarise(avg = mean(expr)) %>%
    pivot_wider(names_from = gene,
                values_from = avg) %>%
    column_to_rownames("disease_class"),

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

  sledai = read_xlsx(path = clinical_file,
                     sheet = "SLEDAI",
                     col_types = c("text",
                                   "text",
                                   "date",
                                   rep("numeric",26),
                                   "date",
                                   rep("numeric",23))) %>%
    clean_names() %>%
    mutate_if(.predicate = is.integer,
              .funs = replace_na,
              0) %>%
    select(-study,
           -form_date) %>%
    left_join(
      read_xlsx(path = "metadata/SLEDAI_BPX.xlsx",
                sheet = "search list",
                .name_repair = janitor::make_clean_names,
                col_types = c("text",
                              "text"))
    ),

  coldata_with_clinical =
    colData(dds_processed) %>%
    as_tibble(rownames = "sample_name") %>%
    left_join(sledai) %>%
    column_to_rownames(var = "sample_name") %>%
    as("DataFrame"),

  # I don't know why and don't care enough to research it, but `colData<-` doesn't seem to return the SummarizedExperiment object
  dds_with_sledai = target({colData(dds_processed) <- coldata_with_clinical; dds_processed}), # he knows this one weird trick...

  dds_with_scores =
    scoreEigengenes(object = dds_with_sledai,
                    module_list = module_list,
                    score_func = 'rsvd') %>%
    scoreEigengenes(object = .,
                    module_list = ldg_modules,
                    score_func = 'rsvd') %>%
    scoreEigengenes(object = .,
                    module_list = metasignature_module,
                    score_func = 'rsvd'),

  annotated_modules = module_annot %>%
    as_tibble(rownames="module") %>%
    filter(type != "Undetermined"),


  ###--- PALETTES ---###
  # we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
  type_pal = paletteer_d("ggsci::category20_d3",
                         n = length(unique(annotated_modules$type))) %>%
    as.character() %>%
    set_names(unique(annotated_modules$type)),

  run_groups =
    colData(dds_with_scores) %>%
    as_tibble() %>%
    pull(run_id) %>%
    unique(),

  run_id_color_set =
    colorRampPalette(
      brewer.pal(12, "Paired"))(length(run_groups)) %>%
    set_names(run_groups),

  chr_pal = c("Y" = "#E41A1C",
              "X" = "#377EB8"),

  sex_pal = c("Male" = "coral3",
              "Female" = "azure2",
              "unk" = "#333333"),

  cluster_pal =
    paletteer::paletteer_d(palette = "ggsci::uniform_startrek",
                           n = length(levels(clusters$cluster))) %>%
    as.character() %>%
    set_names(levels(clusters$cluster)),

  project_pal =
    colorRampPalette(
      brewer.pal(12, "Set1"))(length(levels(annotation_info$project))) %>%
    set_names(levels(annotation_info$project)),

  # disease_class_pal =
  #   paletteer::paletteer_d(palette = "ggthemes::colorblind",
  #                          n = length(unique(annotation_info$disease_class))) %>%
  #   as.character() %>%
  #  set_names(unique(annotation_info$disease_class)),
  disease_class_pal =
    paletteer_d("RColorBrewer::Set1")[seq_along(unique(annotation_info$disease_class))] %>%
    set_names(unique(annotation_info$disease_class)),

  comparison_pal =
    oaColors::oaPalette(length(unique(deg_class$comparison))) %>%
    set_names(unique(deg_class$comparison)),

  group_pal =
    list(
      type_pal,
      chr_pal,
      sex_pal,
      cluster_pal,
      project_pal,
      disease_class_pal,
      comparison_pal) %>%
    set_names(c(
      "type",
      "chr",
      "sex",
      "cluster",
      "project",
      "disease_class",
      "comparison")),

  module_scores =
    colData(dds_with_scores) %>%
    as_tibble(rownames="sample_name") %>%
    select(sample_name,
           matches("^M[[:digit:]]+\\.")) %>%
    column_to_rownames("sample_name"),

  annotated_module_scores =
    colData(dds_with_scores) %>%
    as_tibble(rownames="sample_name") %>%
    select(sample_name,
           one_of(annotated_modules$module)) %>%
    column_to_rownames("sample_name"),

  top_vars = matrixStats::rowSds(vsd_exprs) %>%
    set_names(rownames(vsd_exprs)) %>%
    enframe() %>%
    top_n(20000, value) %>%
    pull(name),

  vsd_top = vsd_exprs[top_vars,] %>% t(),

  sft = pickSoftThreshold(vsd_top,
                          powerVector = c(c(1:10),
                                          seq(from = 12, to=20, by=2)),
                          verbose = 5),

  ###--- WGCNA ---###
  wgcna_modules = blockwiseModules(datExpr = vsd_top,
                                    power = sft$powerEstimate,
                                    maxBlockSize = 20000,
                                    mergeCutHeight = 0.2,
                                    minModuleSize = 20,
                                    pamRespectsDendro = FALSE,
                                    saveTOMs = FALSE,
                                    verbose = 3,
                                    detectCutHeight = 0.995,
                                    TOMDenom = "min",
                                    networkType = "signed hybrid",
                                    reassignThreshold = 1e-6),

  wgcna_module_genes = wgcna_modules$colors %>%
    enframe(name = "gene",
            value = "module"),

  wgcna_hub_genes = chooseTopHubInEachModule(vsd_top,
                                             wgcna_modules$colors,
                                             power = 4,
                                             type = "signed hybrid"),

  wgcna_scores = wgcna_modules$MEs %>%
    as_tibble(rownames="sample_name") %>%
    select(-MEgrey) %>%
    left_join(as_tibble(annotation_info,
                        rownames="sample_name")),

  wgcna_cluster_split = initial_split(data = wgcna_scores %>%
                                        filter(disease_class != "LP") %>%
                                        mutate(disease_class = fct_drop(disease_class)) %>%
                                        select(cluster,
                                               starts_with("ME")),
                                      prop = 0.75,
                                      strata = "cluster"),

  wgcna_cluster_train = training(wgcna_cluster_split),
  wgcna_cluster_test = testing(wgcna_cluster_split),

  wgcna_cluster_rf_cv = train(cluster ~ .,
                       method = "parRF",
                       data = wgcna_cluster_train,
                       trControl = trainControl(method = "repeatedcv",
                                                number = 10,
                                                repeats = 10,
                                                search = "grid",
                                                allowParallel = TRUE),
                       importance=T),

  wgcna_cluster_rf_cv_varImp = varImp(object = wgcna_cluster_rf_cv,
                                      scale = FALSE,
                                      importance=TRUE),

  wgcna_disease_class_split = initial_split(data = wgcna_scores %>%
                                              filter(disease_class != "LP") %>%
                                              mutate(disease_class = fct_drop(disease_class)) %>%
                                              select(disease_class,
                                                     starts_with("ME")),
                                            prop = 0.75,
                                            strata = "disease_class"),

  wgcna_disease_class_train = training(wgcna_disease_class_split),
  wgcna_disease_class_test = testing(wgcna_disease_class_split),

  wgcna_disease_class_rf_cv = train(disease_class ~ .,
                                    method = "parRF",
                                    data = wgcna_disease_class_train,
                                    trControl = trainControl(method = "repeatedcv",
                                                             number = 20,
                                                             repeats = 20,
                                                             search = "grid",
                                                             allowParallel = TRUE)),

  wgcna_disease_class_rf_cv_varImp = varImp(object = wgcna_disease_class_rf_cv,
                                            scale = FALSE,
                                            importance=TRUE),

  ###--- WGCNA-based GSEA ---###
  filtered_wgcna_module_genes =
    wgcna_module_genes %>%
    mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
    filter(!is.na(hugo)),

  MEenriched = future_map_dfr(
    .x = unique(filtered_wgcna_module_genes$module),
    .f = function(i){
      filtered_wgcna_module_genes %>%
        filter(module == i) %>%
        pull(hugo) %>%
        enricher(gene = .,
                 TERM2GENE = c5) %>%
        slot("result") %>%
        mutate(module = i)
      }),

  MEplotting = MEenriched %>%
    filter(p.adjust < 0.05,
           module != "grey") %>%
    mutate(
      GeneRatio = map_dbl(.x = GeneRatio,
                          .f = function(i){
                            j = str_split(i, "/") %>%
                              `[[`(1) %>%
                              as.double()
                            j[[1]]/j[[2]]
                          }),
      ID =  str_replace_all(string = ID, pattern = "_", replacement = " "),
      module = paste0("ME", module)) %>%
    group_by(module) %>%
    top_n(5, GeneRatio) %>%
    sample_n(size = 5,
             replace = TRUE) %>%
    distinct() %>%
    ungroup() %>%
    arrange(module, GeneRatio) %>%
    mutate(order = row_number()),

  viral_transcripts =
    annot %>%
    filter(!str_detect(string = transcript,
                       pattern = "^ENST")) %>%
    pull(gene_name) %>%
    intersect(rownames(vsd_exprs)),

  detected_viral_transcripts =
    counts(dds_with_scores) %>%
    t() %>%
    as_tibble(rownames = "sample_name") %>%
    select(sample_name,
           one_of(viral_transcripts)) %>%
    pivot_longer(-sample_name,
                 names_to = "transcript",
                 values_to = "counts") %>%
    group_by(transcript) %>%
    summarise(total_counts = sum(counts)) %>%
    filter(total_counts > 0) %>%
    pull(transcript),

  viral_exprs =
    vsd_exprs[detected_viral_transcripts,] %>%
    t() %>%
    as_tibble(rownames="sample_name"),

  ifn_modules =
    annotated_modules %>%
    filter(type == "Interferon") %>%
    pull(module),

  inflame_modules =
    annotated_modules %>%
    filter(type == "Inflammation") %>%
    pull(module),

  module_scores_with_viral =
    colData(dds_with_scores) %>%
    as_tibble(rownames='sample_name') %>%
    select(sample_name,
           disease_class,
           matches("^M[[:digit:]]+"),
           one_of(names(ldg_modules)),
           mg) %>%
    inner_join(wgcna_scores) %>%
    inner_join(viral_exprs) %>%
    inner_join(clusters),

  report = rmarkdown::render(
    input = knitr_in("markdown/report.rmd"),
    output_file = file_out("markdown/report.html"),
    output_dir = "markdown",
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

  wgcna_modules$MEs %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/wgcna_eigengene_scores.csv")),

  wgcna_module_genes %>%
    write_csv(path = file_out("data/wgcna_genes.csv")),

  colData(dds_with_scores) %>%
    as_tibble(rownames="sample_name") %>%
    select(sample_name,
           matches("^M[[:digit:]]+"),
           one_of(names(ldg_modules))) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("data/module_scores.csv")),

  disp_plot = plot_dispersion_estimate(dds_with_scores),

  save(annotation_info,
       final_md,
       group_pal,
       module_scores_with_viral,
       ISGs,
       pca_results,
       degs,
       umap_results,
       vsd_exprs,
       wgcna_modules,
       file = file_out("data/for_shiny_vis.RData"))
)
