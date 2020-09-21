#For ABC-MESS
analysis_plan = drake_plan(

  study_metadata = read_excel(path = metadata_file,
                              sheet = "main",
                              col_types = c("numeric",
                                            rep("text",10),
                                            "numeric",
                                            rep("text",6),
                                            rep("numeric", 3),
                                            "logical",
                                            "logical"),
                              trim_ws = TRUE,
                              na = "n/a") %>%
    clean_names() %>%
    mutate(sample_name = janitor::make_clean_names(sample_name, case = "all_caps"),
           disease_class = str_remove(
             string = recode(disease_class, "Unaffected Control" = "Control"),
             pattern = " Patient"),
           mess = replace_na(data = mess, replace = FALSE),
           abc_or_mess_or_control = replace_na(data = abc_or_mess_or_control, replace = FALSE)
    ) %>%
    filter(initial_concentration_ng_ul > initial_concentration_threshold),

  non_project_controls =
    study_metadata %>%
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
           rin,
           subject_ref,
           visit_ref),

  abc_mess_samples =
    study_metadata %>%
    filter(
      abc_or_mess_or_control == TRUE,
      project %nin% projects_to_exclude
    ) %>%
    # Select the portions of the metadata that are useful:
    select(
      sample_name,
      disease_class,
      project,
      run_id,
      sample_alias,
      sex,
      race_code,
      initial_concentration_ng_ul,
      final_concentration_ng_ul,
      rin,
      subject_ref,
      visit_ref
      ),

  medication_md =
    read_excel(
      path = "metadata/abcmess.xlsx",
      sheet = "ABCmed",
      col_types = c(
        rep("text",5),
        "date",
        "numeric",
        "text",
        rep("numeric", 5),
        rep("text", 2),
        "logical",
        rep("text",5)
      ),
      .name_repair = janitor::make_clean_names) %>%
    set_names(str_remove_all(string = names(.), pattern = "_mg_[[:graph:]]+")) %>%
    mutate(milestone =
             tolower(milestone) %>%
             str_replace(
               pattern = " ",
               replacement = "_"
             ) %>%
             as_factor(),
           sledai_group =
             tolower(sledai_group) %>%
             as_factor(),
           abatacept =
             replace_na(
               data = abatacept,
               replace = FALSE),
           treatment_status =
             tolower(treatment_status) %>%
             as_factor(),
           depomedrol_methylprednisone =
             recode(depomedrol_methlyprednisone,
                    "40 not for SLE" = "40",
                    "160 10/24/16" = "160",
                    "depo dose pack 6/22-6/27/17" = "NA",
                    "160 7/19/17" = "160",
                    "24-4mg medrol dose pack 9/2/17-9/7/17" = "24",
                    "160  (540mg  4/19/18 to 4/26/18)" = "160",
                    "120 7/2/18" = "120") %>%
             as.integer(),
           prednisone =
             recode(prednisone,
                    "5-7.5 EOD" = "5",
                    "5  3/7/16" = "5",
                    "30-10mg 9/1/16-9/14/16" = "10"
             ) %>%
             as.integer(),
           other_medications =
             replace_na(
               data = other_medications,
               replace = "none"
             ) %>%
             as_factor(),
           responders =
             tolower(responders) %>%
             as_factor()
    ) %>%
    select(-id_2, -treatment_status_2, -id_3, -note, -depomedrol_methlyprednisone),

  md =
    study_metadata %>%
    filter(project %in% (projects_to_include %||% unique(.data$project)),
           project %nin% projects_to_exclude) %>%
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
           rin,
           subject_ref,
           visit_ref) %>%
    bind_rows(non_project_controls) %>%
    bind_rows(abc_mess_samples) %>%
    rename(ethnicity = race_code) %>%
    mutate(
      project = as_factor(project),
      disease_class = as_factor(disease_class),
      run_id = factor(run_id),
      initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200),
      disease_class = as_factor(disease_class),
      sex = as_factor(sex),
      ethnicity = as_factor(ethnicity)) %>%
    filter(
      disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
      disease_class %nin% disease_classes_to_exclude) %>%
    distinct() %>%
    left_join(medication_md),

  tx_sample_names = dir(path = seq_file_directory,
                        pattern = "quant.sf.gz",
                        recursive = TRUE,
                        full.name = TRUE) %>%
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>%
    str_split(pattern = "/") %>%
    map_chr(~pluck(.x, length(.x)-1)) %>%
    str_remove(pattern = '(_[L|S][[:digit:]]+)+') %>%
    janitor::make_clean_names(case = "all_caps"),

  tx_files = dir(path = seq_file_directory,
                 pattern = "quant.sf.gz",
                 recursive = TRUE,
                 full.name = TRUE) %>%
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>%
    set_names(tx_sample_names) %>%
    purrr::discard(.p = is.na(match(names(.), md$sample_name))),

  final_md = md %>%
    filter(
      sample_name %in% names(tx_files),
      str_detect(
        string = sample_name,
        pattern = "_2$",
        negate = TRUE
      )
    ) %>%
    mutate(
      disease_class = fct_relevel(disease_class, {{control_group}}),
      run_id = fct_drop(run_id)) %>%
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

  corrected_counts = ComBat_seq(counts = counts(dds_import),
                                batch = fct_drop(colData(dds_import)[[batch_variable]]),
                                group = colData(dds_import)[[comparison_grouping_variable]]) %>%
    `storage.mode<-`("integer"),

  dds_import_combat =
    DESeqDataSetFromMatrix(
      countData = corrected_counts,
      colData = colData(dds_import),
      design = study_design
    ),

  dds_filtered = dds_import_combat %>%
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

  dds =
    DESeq(dds_qc,
          parallel = TRUE),
  sva_res = calc_sva(dds = dds, model_design = "disease_class", n.sva = num_sva),
  dds_processed = sva_res$dds,
  sva_graph_data = sva_res$sva,
  vsd = vst(dds_processed),
  vsd_exprs = assay(vsd),

  sample_dists = vsd_exprs %>% t() %>% parallelDist::parallelDist(),
  #fig.width=12, fig.height=9
  sampleDistMatrix = as.matrix(sample_dists),

  # Minus the noise, actual differences and similarities are now apparent.
  sample_dendrogram = sample_dists %>% hclust() %>% as.dendrogram(),

  sample_cluster_info = ident_clusters(annotated_module_scores, K.max = 20),

  clusters = sample_cluster_info$clusters %>% mutate(cluster = as_factor(cluster)),

  annotation_info = as.data.frame(colData(dds_with_scores))[,c("disease_class",
                                                               "sex")] %>%
    as_tibble(rownames = "sample_name") %>%
    left_join(clusters) %>%
    column_to_rownames(var="sample_name"),

  # Variation can also be examined in reduced dimensional space by PCA or UMAP:
  pca_results =
    irlba::prcomp_irlba(x = vsd_exprs) %>%
    pluck("rotation") %>%
    as_tibble() %>%
    mutate(sample_name = colnames(vsd_exprs)) %>%
    inner_join(as_tibble(colData(dds_processed),
                         rownames="sample_name")) %>%
    inner_join(clusters),

  umap_results =
    umap(
      t(vsd_exprs),
      n_threads = detectCores(),
      n_sgd_threads = detectCores(),
      verbose = TRUE,
      n_components = 3
    ) %>%
    as_tibble(.name_repair = "unique") %>%
    set_names(c("umap_1",
                "umap_2",
                "umap_3")) %>%
    mutate(sample_name = colnames(vsd_exprs)) %>%
    inner_join(as_tibble(colData(dds_processed),
                         rownames="sample_name")) %>%
    inner_join(clusters),

  comparison_results_list =
    resultsNames(dds_processed) %>%
    keep(
      str_detect(
        string = .,
        pattern = comparison_grouping_variable)),
  res = map(comparison_results_list, function(i){
    lfcShrink(dds_processed,
              coef = i,
              parallel = TRUE,
              type = "apeglm")
  }) %>%
    set_names(map_chr(comparison_results_list, str_remove, pattern = paste0(comparison_grouping_variable, "_"))),

  down_tables = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        !is.na(padj) & padj <= 0.05,
        log2FoldChange < 0
      ) %>%
      mutate(
        log2FoldChange = -log2FoldChange
      ) %>%
      mutate_at(
        .vars = vars(-gene),
        .funs = list(~signif(x = ., digits =  2))
      ) %>%
      top_n(
        n = 25,
        wt = log2FoldChange
      ) %>%
      arrange(
        desc(
          log2FoldChange
        )
      )
  }) %>%
    set_names(
      nm = map_chr(
        .x = comparison_results_list,
        .f = str_remove,
        pattern = str_glue("{comparison_grouping_variable}_")
      )
    ) %>%
    keep(~ nrow(.x) > 0),

  up_tables = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        !is.na(padj) & padj <= 0.05,
        log2FoldChange > 0
      ) %>%
      mutate_at(
        vars(-gene),
        list(~signif(., 2)
        )
      ) %>%
      top_n(
        n = 25,
        wt = log2FoldChange
      ) %>%
      arrange(
        desc(
          log2FoldChange
        )
      )
  }) %>%
    set_names(
      nm = map_chr(
        .x = comparison_results_list,
        .f = str_remove,
        pattern = str_glue("{comparison_grouping_variable}_")
      )
    ) %>%
    keep(~ nrow(.x) > 0),

  degs = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        padj < 0.05
      ) %>%
      filter(
        abs(
          log2FoldChange
        ) >= 0.5
      ) %>%
      pull(gene)
  }) %>%
    set_names(
      map_chr(
        .x = comparison_results_list,
        .f = str_remove,
        pattern = str_glue("{comparison_grouping_variable}_")
      )
    ),

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
                select(name, sex, ethnicity, disease_class)) %>%
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

  banchereau_gene_module_tbl =
    module_list %>%
    enframe(name = "module",
            value = "gene") %>%
    unnest(cols = "gene"),

  ldg_gene_module_tbl =
    ldg_modules %>%
    enframe(name = "module",
            value = "gene") %>%
    unnest(cols = "gene"),

  mg_gene_module_tbl =
    metasignature_module %>%
    enframe(name = "module",
            value = "gene") %>%
    unnest(cols = "gene"),

  dds_with_scores =
    scoreEigengenes(object = dds_processed,
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

  cluster_pal = ifelse(
    test = length(levels(clusters$cluster)) > 12,
    yes = list(
      colorRampPalette(
        paletteer_d(
          palette = "ggthemes::calc",
          n = 12
        )
      )(
        length(
          levels(
            clusters$cluster
          )
        )
      )
    ),
    no = list(
      paletteer_d(
        palette = "ggthemes::calc",
        n = length(
          levels(
            clusters$cluster
          )
        )
      )
    )
  ) %>%
    unlist() %>%
    as.character() %>%
    set_names(levels(clusters$cluster)),

  project_pal =
    colorRampPalette(
      brewer.pal(12, "Set1"))(length(levels(annotation_info$project))) %>%
    set_names(levels(annotation_info$project)),

  number_disease_classes = length(unique(annotation_info$disease_class)),

  disease_class_pal =
    if_else(
      number_disease_classes > 2,
      list(RColorBrewer::brewer.pal(number_disease_classes, "Set1")),
      list(c("black", "grey75"))
    ) %>%
    unlist() %>%
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
  wgcna_modules =
    blockwiseModules(
      datExpr = vsd_top,
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
      reassignThreshold = 1e-6
    ),

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
                                        mutate(disease_class = fct_drop(disease_class)) %>%
                                        select(cluster,
                                               starts_with("ME")),
                                      prop = 0.75,
                                      strata = "cluster"),

  wgcna_cluster_train = training(wgcna_cluster_split),
  wgcna_cluster_test = testing(wgcna_cluster_split),

  wgcna_cluster_rf_cv =
    train(
      cluster ~ .,
      method = "parRF",
      data = wgcna_cluster_train,
      trControl =
        trainControl(
          method = "repeatedcv",
          number = 10,
          repeats = 10,
          search = "grid",
          allowParallel = TRUE
        ),
      importance=T
    ),

  wgcna_cluster_rf_cv_varImp = varImp(object = wgcna_cluster_rf_cv,
                                      scale = FALSE,
                                      importance=TRUE),

  wgcna_disease_class_split = initial_split(data = wgcna_scores %>%
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

  ###---module_classification---###
  module_scores_with_md = module_scores %>%
    as_tibble(rownames="sample_name") %>%
    left_join(as_tibble(annotation_info,
                        rownames="sample_name")),

  module_cluster_split = initial_split(data = module_scores_with_md %>%
                                         filter(disease_class != "LP") %>%
                                         mutate(disease_class = fct_drop(disease_class)) %>%
                                         select(cluster,
                                                matches("^M[0-9]+")),
                                       prop = 0.75,
                                       strata = "cluster"),

  module_cluster_train = training(module_cluster_split),
  module_cluster_test = testing(module_cluster_split),

  module_cluster_rf_cv = train(cluster ~ .,
                               method = "parRF",
                               data = module_cluster_train,
                               trControl = trainControl(method = "repeatedcv",
                                                        number = 10,
                                                        repeats = 10,
                                                        search = "grid",
                                                        allowParallel = TRUE),
                               importance=T),

  module_cluster_rf_cv_varImp = varImp(object = module_cluster_rf_cv,
                                       scale = FALSE,
                                       importance=TRUE),

  module_disease_class_split = initial_split(data = module_scores_with_md %>%
                                               filter(disease_class != "LP") %>%
                                               mutate(disease_class = fct_drop(disease_class)) %>%
                                               select(disease_class,
                                                      matches("^M[0-9]+")),
                                             prop = 0.75,
                                             strata = "disease_class"),

  module_disease_class_train = training(module_disease_class_split),
  module_disease_class_test = testing(module_disease_class_split),

  module_disease_class_rf_cv = train(disease_class ~ .,
                                     method = "parRF",
                                     data = module_disease_class_train,
                                     trControl = trainControl(method = "repeatedcv",
                                                              number = 10,
                                                              repeats = 10,
                                                              search = "grid",
                                                              allowParallel = TRUE),
                                     importance=T),

  module_disease_class_rf_cv_varImp = varImp(object = module_disease_class_rf_cv,
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
                 TERM2GENE = c5)%>%
        {if(!is.null(.)) mutate(slot(.,'result'), module = i) } #if_else does not work because it tries to evaulate both the true and false results, which doesn't work with `slot(NULL, "results")`
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

  cd72_and_steroids = final_md %>%
    filter(sample_name %in% colnames(vsd_exprs)) %>%
    left_join(t(vsd_exprs)[,"CD72", drop = FALSE] %>%
                as_tibble(rownames = "sample_name")) %>%
    select(prednisone, depomedrol_methylprednisone, disease_class, CD72, sample_name) %>%
    mutate(prednisone = replace_na(prednisone, 0),
           depomedrol_methylprednisone = replace_na(depomedrol_methylprednisone, 0)),

  cd72_prednisone_correlation =
    cd72_and_steroids %>%
    group_by(disease_class) %>%
    cor_test(vars = CD72, vars2 = prednisone),

  cd72_depomedrol_correlation =
    cd72_and_steroids %>%
    group_by(disease_class) %>%
    cor_test(vars = CD72, vars2 = depomedrol_methylprednisone),


  report = rmarkdown::render(
    input = knitr_in("markdown/report.rmd"),
    output_file = file_out("results/report.html"),
    output_dir = "results",
    quiet = TRUE),

  qc_report = rmarkdown::render(
    input = knitr_in("markdown/qc_report.rmd"),
    output_file = file_out("results/qc_report.html"),
    output_dir = "results",
    quiet = TRUE),

  supplemental_report = rmarkdown::render(
    input = knitr_in("markdown/supplemental_report.rmd"),
    output_file = file_out("results/supplemental_report.html"),
    output_dir = "results",
    quiet = TRUE),

  vsd_exprs %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("results/filtered_normalized_stabilized_expression.csv")),
  colData(dds_with_scores) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("results/filtered_metadata.csv")),
  counts(dds_with_scores) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("results/filtered_transcript_counts.csv")),
  res %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("results/log_fold_changes.csv")),

  wgcna_modules$MEs %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("results/wgcna_eigengene_scores.csv")),

  wgcna_module_genes %>%
    write_csv(path = file_out("results/wgcna_genes.csv")),

  colData(dds_with_scores) %>%
    as_tibble(rownames="sample_name") %>%
    select(sample_name,
           matches("^M[[:digit:]]+"),
           one_of(names(ldg_modules))) %>%
    as.data.frame() %>%
    data.table::fwrite(row.names = TRUE,
                       file = file_out("results/module_scores.csv")),

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
       file = file_out("results/for_shiny_vis.RData"))
)
