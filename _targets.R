library(targets)
library(tarchetypes)
source("code/packages.R")
source("code/plan/01_import_funcs.R")
source("code/plan/02_filtering_funcs.R")
source("code/plan/04_module_funcs.R")
source("code/plan/06_dimensional_reduction_funcs.R")
source("code/plan/07_differential_expression_funcs.R")
source("code/plan/10_viral_transcript_funcs.R")
source("code/plan/11_stats_testing_funcs.R")
source("code/plan/98_palettes_funcs.R")

source("project_options.R")
options(tidyverse.quiet = TRUE)

list(
  tar_target(
    md,
    import_metadata(
      metadata_file = metadata_file,
      metadata_sheet = "main",
      extra_controls_metadata_file = main_sample_list,
      extra_controls_metadata_sheet = "main",
      groups_to_include = project_groups_to_include,
      groups_to_exclude = project_groups_to_exclude)
    ),
  tar_target(
    tx_files,
    import_counts(
      seq_file_directory = seq_file_directory
    )
  ),
  tar_target(
    dds,
    create_initial_deseq_dataset(
      metadata = md,
      tx_files = tx_files,
      study_design = study_design,
      comparison_group = comparison_grouping_variable,
      control_group = control_group
    )
  ),

  tar_target(
    dds_filtered,
    filter_counts(
      dds = dds,
      min_counts = 1, 
      removal_pattern = "^RNA5"
      )
    ),

  tar_target(
    outlier_qc
    remove_outliers(
      dds = dds_filtered,
      pc1_zscore_cutoff = pc1_zscore_threshold,
      pc2_zscore_cutoff = pc2_zscore_threshold
      )
  ),

  tar_target(
    dds_qc,
    DESeq(
      outlier_qc$dds,
      parallel = TRUE,
      BPPARAM = BPPARAM
    )
  ),
  
  tar_target(
    sva_res,
    calc_sva(
      dds = dds,
      model_design = comparison_grouping_variable,
      n.sva = num_sva
      )
  ),

  tar_target(
    vsd,
    vst(sva_res$dds)
  ),

  tar_target(
    vsd_exprs,
    assay(vsd)
  ),

  tar_target(
    dds_with_scores,
    scoreEigengenes(
      object = dds_processed,
      module_list = banchereau_modules,
      score_func = 'rsvd'
      ) %>%
    scoreEigengenes(
      object = .,
      module_list = ldg_modules,
      score_func = 'rsvd'
      ) %>%
    scoreEigengenes(
      object = .,
      module_list = metasignature_module,
      score_func = 'rsvd'
      )
  ),

  tar_target(
    module_scores,
    extract_module_scores(
      dds = dds_with_scores, 
        names(banchereau_modules),
        names(ldg_modules),
        names(metasignature_module)
        )
  ),

  tar_target(
    disp_plot,
    plot_dispersion_estimate(dds_with_scores)
  ) ,

  tar_target(
    annotated_modules,
    filter(
      .data = module_annotation,
      type != "Undetermined"
      )
  ),

  tar_target(
    annotated_mod_list,
    mutate(
      .data = annotated_modules,
      module_type =
        paste(
          module,
          type,
          sep = " - "
        )
    ) %>%
    select(-type) %>%
    deframe()
  ),

  tar_target(
    annotated_module_scores,
    select(
      .data = module_scores,
      sample_name,
      one_of(annotated_modules$module)
      )
  ),

  tar_target(
    module_tbl,
    create_module_table(
      banchereau_modules,
      ldg_modules,
      metasignature_module
    )
  ),

  tar_target(
    module_ISGs,
    extract_module_genes(
      module_table      = module_tbl,
      module_annotation = "Interferon"
    )
  ),
  
  tar_target(
    sample_dists
    parallelDist(t(vsd_exprs))
  ),

  tar_target(
    sample_dendrogram,
    as.dendrogram(hclust(sample_dists))
  ),

  tar_target(
    sample_cluster_info,
    ident_clusters(
      column_to_rownames(
        annotated_module_scores,
        "sample_name"
        ),
      K.max = 20
    )
  ),

  tar_target(
    clusters,
    mutate(
      .data = sample_cluster_info$clusters,
      cluster = as_factor(cluster)
    )
  ),

  tar_target(
    study_md,
    left_join(
      as_tibble(
        colData(dds_with_scores),
        rownames = "sample_name"),
      clusters
    )
  ),

  tar_target(
    annotation_info,
    select(.data = study_md,
          disease_class,
          sex,
          cluster,
          sample_name) %>%
      column_to_rownames(var = "sample_name")
  ),

  tar_target(
    pca_results,
    run_pca(
      expr_data = vsd_exprs,
      metadata = 
        as_tibble(
          x = colData(dds_with_scores),
          rownames="sample_name"
          ),
      cluster_info = clusters
    )
  ),

  tar_target(
    umap_results,
    run_umap(
      expr_data = vsd_exprs,
      metadata = 
        as_tibble(
          x = colData(dds_with_scores),
          rownames = "sample_name"
      ),
      cluster_info = clusters
    )
  ),

  tar_target(
    comparison_results_list,
    keep(
      resultsNames(dds_with_scores),
      str_detect(
        string = resultsNames(dds_with_scores),
        pattern = comparison_grouping_variable
      )
    )
  ),

  tar_target(
    name = res,
    command = 
      create_results_list(
        comparison_list = comparison_results_list,
        dds = dds,
        comparison_grouping_variable = comparison_grouping_variable
      )
  ),

  tar_target(
    name = down_tables,
    command = create_downregulation_tables(
      results = res,
      comparison_list = comparison_results_list,
      grouping_variable = comparison_grouping_variable
    )
  ),

  tar_target(
    name = up_tables,
    command = create_downregulation_tables(
      results = res,
      comparison_list = comparison_results_list,
      grouping_variable = comparison_grouping_variable
    )
  ),

  tar_target(
    name = degs,
    command = extract_de_genes(
      results = res,
      comparison_list = comparison_results_list,
      grouping_variable = comparison_grouping_variable
    )
  ),

  tar_target(
    name = deg_class,
    command = group_degs(degs)
  )

  tar_target(
    name = deg_means,
    command = calc_deg_means(
      exprs = vsd_exprs,
      deg_class = deg_class,
      metadata = as_tibble(
        colData(dds_with_scores),
        rownames = "name"
      ) %>%
        select(
          name,
          sex,
          ethnicity,
          disease_class
        ),
      grouping_variable = disease_class
    )
  ),

  tar_target(
    name = ISGs,
    command = intersect(
      c(
        "STAT1", "ADAR", "ABCE1", "RNASEL", "TYK2", "IFNAR1",
        "IFNB1", "STAT2", "IFNAR2", "JAK1", "SAMHD1",
        "SOCS1", "SOCS3", "STAT1", "ISG20", "IFITM3",
        "IFITM1", "IRF9", "ISG15", "IFI6", "IFIT3", "USP18",
        "IP6K2", "PSMB8", "IFIT1", "IRF4", "IRF5", "IRF1", "IRF3",
        "IRF6", "IRF8", "IRF2", "IRF7", "IFITM2", "XAF1", "IFI27",
        "GBP2", "RSAD2", "MX2", "MX1", "IFIT2", "IFI35", "BST2",
        "OAS1", "OASL", "OAS2", "OAS3", "PTPN6", "PTPN11"
        ),
      rownames(vsd_exprs)
      )
  ),
  
  tar_target(
    name = vsd_top,
    command = top_variable_genes(
      exprs = vsd_exprs,
      n = 20000
    )
  )

  tar_target(
    name = sft,
    command = pickSoftThreshold(
        data = vsd_top,
        powerVector =
          c(
            seq(10),
            seq(
              from = 12,
              to = 30,
              by = 1)
          ),
        verbose = 5
      )
  )

  tar_target(
    name = wgcna_modules,
    command = blockwiseModules(
      datExpr =vsd_top,
      power = find_softPower(sft),
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
    )
  )

  tar_target(
    name = wgcna_module_genes,
    command =
      enframe(
        x = wgcna_modules$colors
        name = "gene",
        value = "module"
      )
  )

  tar_target(
    name = wgcna_modules,
    command = unique(pull(wgcna_module_genes, module))
  )

  tar_target(
    name = wgcna_hub_genes,
    command = 
      chooseTopHubInEachModule(
        datExpr = vsd_top_float,
        colorh = wgcna_modules$colors,
        power = 4,
        type = "signed hybrid"
      )
  )

  tar_target(
    name = wgcna_scores, 
    command = 
      left_join(
        x = select(.data = as_tibble(x = wgcna_modules$MEs, rownames="sample_name"), -MEgrey),
        y = as_tibble(x = annotation_info, rownames="sample_name")
        )
  )

  tar_target(
    name = wgcna_cluster_split,
    command = 
      initial_split(
        data = wgcna_scores %>%
          mutate(disease_class = fct_drop(disease_class)) %>%
          select(
            cluster,
            starts_with("ME")
          ),
        prop = 0.75,
        strata = "cluster"
      )
  )

  tar_target(
    name = wgcna_cluster_train,
    command =  training(wgcna_cluster_split)
  )

  tar_target(
    name = wgcna_cluster_test,
    command = testing(wgcna_cluster_split)
  )

  tar_target(
    name = wgcna_cluster_rf_cv,
    command =
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
        importance=TRUE
      )
  )

  tar_target(  
    name = wgcna_cluster_rf_cv_varImp =
    command =
      varImp(
        object = wgcna_cluster_rf_cv,
        scale = FALSE,
        importance = TRUE
      ))

  tar_target(
    name = wgcna_disease_class_split,
    command =
      initial_split(
        data = 
          select(
            .data = mutate(
              .data = wgcna_scores,
              disease_class = fct_drop(disease_class)
              ),
            disease_class,
            starts_with("ME")
          ),
        prop = 0.75,
        strata = "disease_class"
      )
  )

  tar_target(
    name = wgcna_disease_class_train,
    command = training(wgcna_disease_class_split)
  )

  tar_target(
    name = wgcna_disease_class_test,
    command = testing(wgcna_disease_class_split)
  ),

  tar_target(
    name = wgcna_disease_class_rf_cv,
    commans =
      train(
        form = disease_class ~ .,
        method = "parRF",
        data = wgcna_disease_class_train,
        trControl =
          trainControl(
            method = "repeatedcv",
            number = 20,
            repeats = 20,
            search = "grid",
            allowParallel = TRUE
          )
      )
  ),

  tar_target(
    name = wgcna_disease_class_rf_cv_varImp,
    command = 
      varImp(
        object = wgcna_disease_class_rf_cv,
        scale = FALSE,
        importance=TRUE
      )
  )

  tar_target(
    name = filtered_wgcna_module_genes,
    command = 
      filter(
        .data =
          mutate(
            .data = wgcna_module_genes,
            hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]
          ),
        !is.na(hugo)
      )
  ),

  tar_target(
    name = MEenriched,
    command = module_gsea(
      module_genes = filtered_wgcna_module_genes,
      module_of_interest = i
    ),
    pattern = map(wgcna_modules)
  ),

  tar_target(
    name = MEplotting,
    command = module_gsea_plots(enriched_genes = MEenriched)
  ),

  tar_target(
    name = module_scores_with_md,
    command = 
      left_join(
        x = module_scores,
        y = as_tibble(
          x = annotation_info,
          rownames="sample_name"
          )
      )
  ),

  tar_target(
    name = module_cluster_split,
    command =
      initial_split(
        data =
          select(
            .data = 
              mutate(
                .data = module_scores_with_md,
                disease_class = fct_drop(disease_class)
              ),
            cluster,
            one_of(names(banchereau_modules))
          ),
        prop = 0.75,
        strata = "cluster"
      )
  ),

  tar_target(
    name = module_cluster_train,
    command = training(module_cluster_split)
  ),

  tar_target(
    name = module_cluster_test,
    command = testing(module_cluster_split)
  ),

  tar_target(
    name = module_cluster_rf_cv,
    command = 
      train(
        form = cluster ~ .,
        method = "parRF",
        data = module_cluster_train,
        trControl =
          trainControl(
            method = "repeatedcv",
            number = 10,
            repeats = 10,
            search = "grid",
            allowParallel = TRUE
          ),
        importance = TRUE
      )
    ),

  tar_target(
    name = module_cluster_rf_cv_varImp,
    command = 
      varImp(
        object = module_cluster_rf_cv,
        scale = FALSE,
        importance = TRUE
      )
  ),

  tar_target(
    name = module_disease_class_split,
    command =
      initial_split(
        data =
          select(
            .data = 
              mutate(
                .data = module_scores_with_md,
                disease_class = fct_drop(disease_class)
                ),
            disease_class,
            one_of(names(banchereau_modules))
          ),
        prop = 0.75,
        strata = "disease_class"
      )
  ),

  tar_target(
    name = module_disease_class_train
    command = training(module_disease_class_split)
  ),

  tar_target(
    name = module_disease_class_test,
    command = testing(module_disease_class_split)
  ),

  tar_target(
    name = module_disease_class_rf_cv,
    command =
      train(
        form = disease_class ~ .,
        method = "parRF",
        data = module_disease_class_train,
        trControl =
          trainControl(
            method = "repeatedcv",
            number = 10,
            repeats = 10,
            search = "grid",
            allowParallel = TRUE
          ),
        importance = TRUE
      )
  ),

  tar_target(
    name = module_disease_class_rf_cv_varImp,
    command =
      varImp(
        object = module_disease_class_rf_cv,
        scale = FALSE,
        importance = TRUE
      )
  ),

  tar_target(
    name = annot,
    command = load_gene_annotations(),
    format = "file"
  ),

  tar_target(
    name = viral_exprs,
    command =
      extract_viral_expression(
        annotations = annot,
        exprs = vsd_exprs,
        dds = dds_with_scores
        )
  ),

  tar_target(
    name = ifn_modules,
    command = 
      pull(
        .data =
          filter(
            .data = annotated_modules,
            type == "Interferon"
            ),
        var = module
      )
  )

  tar_target(
    name = inflame_modules,
    command = 
      pull(
        .data =
          filter(
            .data = annotated_modules,
            type == "Inflammation"
            ),
        var = module
      )
  ),

  tar_target(
    name = module_scores_with_viral =
    command =
      select(
        .data = study_md,
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
  ),

  tar_target(
    name = annotated_module_scores_with_cluster_class,
    command = 
      select(
        .data = 
        mutate(
          .data =module_scores_with_viral,
          cluster = as_factor(cluster)
        ),
        cluster,
        disease_class,
        one_of(annotated_modules$module)
      )
  ),

  tar_target(
    name = renamed_annotated_module_scores,
    command = 
      set_names(
        nm = 
          drake_recode(
            target_list = names(annotated_module_scores_with_cluster_class)
            thing_to_unquote_splice = annotated_mod_list),
        x  = annotated_module_scores_with_cluster_class
      )
  ),

  tar_target(
    name = annotated_module_scores_pivot,
    command = 
      pivot_longer(
        data      = renamed_annotated_module_scores,
        cols      = starts_with("M"),
        names_to  = "module",
        values_to = "score"
    )
  ),

  tar_target(
    name = annotated_module_stats_by_cluster,
    command =
      modules_compare_with_stats(
        module_score_table = annotated_module_scores_pivot,
        comparison = cluster
      )
  ),

  tar_target(
    name = annotated_module_stats_by_disease,
    command =
      modules_compare_with_stats(
        module_score_table = annotated_module_scores_pivot,
        comparison = disease_class
      )
  ),

  tar_target(
    name = module_scores_pivot,
    command = pivot_module_scores(module_scores = module_scores_with_viral)
  ),

   tar_target(
    name = module_stats_by_cluster,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        comparison = cluster
      )
  ),

   tar_target(
    name = module_stats_by_disease,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        comparison = disease_class
      )
  ),

  tar_target(
    name = module_scores_with_viral_by_cluster,
    command = 
    mutate(
      .data = 
      pivot_longer(
        data =
          select(
            .data = module_scores_with_viral,
            cluster,
            matches("^ME")
          ),
        -cluster,
        names_to     = "module",
        values_to    = "score"
        ),
      cluster = as_factor(cluster)
    )
  ),

  tar_target(
    name = module_scores_with_viral_by_cluster_stats,
    command = 
      modules_compare_with_stats(
        module_score_table = module_scores_with_viral_by_cluster,
        comparison         = cluster
      )
  ),

  tar_target(
    name = module_scores_with_viral_by_disease,
    command = 
    mutate(
      .data = 
      pivot_longer(
        data =
          select(
            .data = module_scores_with_viral,
            disease_class,
            matches("^ME")
          ),
        -cluster,
        names_to     = "module",
        values_to    = "score"
        ),
      disease_class = as_factor(disease_class)
    )
  ),

  tar_target(
    name = module_scores_with_viral_by_disease_stats,
    command = 
      modules_compare_with_stats(
        module_score_table = module_scores_with_viral_by_disease,
        comparison         = disease_class
      )
  ),

  tar_target(
    name = group_pal
    command = 
      create_palettes(
        annotated_modules = annotated_modules,
        clusters = clusters,
        annotation_info = annotation_info,
        deg_class = deg_class
      )
  ),

  tar_render(
    name          = primary_report,
    path          = "analysis/report.rmd",
    params        = list(
      set_title   = "Initial COVID OSCTR samples RNAseq Analysis",
      set_author  = "Miles Smith"
      ),
    ),

  tar_render(
    name          = qc_report,
    path          =  "analysis/qc_report.rmd",
    params        = list(
      set_title   = "Initial COVID OSCTR samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
  )
)