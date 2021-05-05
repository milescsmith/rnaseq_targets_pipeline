library(targets)
library(tarchetypes)

source("code/project_parameters.R")

source("code/plan/01_import_funcs.R")
source("code/plan/02_filtering_funcs.R")
source("code/plan/04_module_funcs.R")
source("code/plan/06_dimensional_reduction_funcs.R")
source("code/plan/07_differential_expression_funcs.R")
source("code/plan/08_WGCNA_funcs.R")
source("code/plan/10_viral_transcript_funcs.R")
source("code/plan/11_stats_testing_funcs.R")
source("code/plan/13_pathways.R")
source("code/plan/98_palettes_funcs.R")
source("code/plan/99_output_funcs.R")

options(tidyverse.quiet = TRUE)

tar_option_set(
  packages = c(
    "caret",
    "cluster",
    "clusterProfiler",
    "corrplot",
    "cowplot",
    "data.table",
    "DESeq2",
    "drake",
    "factoextra",
    "flextable",
    "formattable",
    "furrr",
    "ggbeeswarm",
    "ggforce",
    "ggplotify",
    "ggpubr",
    "ggradar",
    "ggrepel",
    "ggtext",
    "gtools",
    "here",
    "HGNChelper",
    "inspectdf",
    "irlba",
    "janitor",
    "kableExtra",
    "knitr",
    "magrittr",
    "matrixStats",
    "moduleScoreR",
    "oaColors",
    "paletteer",
    "pheatmap",
    "plotly",
    "randomForest",
    "RColorBrewer",
    "readxl",
    "rlang",
    "rmarkdown",
    "rstatix",
    "scales",
    "stats",
    "sva",
    "tidyHeatmap",
    "tidymodels",
    "tidyverse",
    "tximport",
    "uwot",
    "viridis",
    "WGCNA"
  )
)

list(
  tar_target(
    name       = raw_metadata,
    command    = project_params[["metadata_file"]],
    format     = "file",
    deployment = "main"
  ),

  tar_target(
    name       = md,
    command    =
      import_metadata(
        metadata_file = raw_metadata,
        projects_to_include     = project_params[["projects_to_include"]],
        projects_to_exclude     = project_params[["projects_to_exclude"]],
        study_groups_to_include = project_params[["study_groups_to_include"]],
        study_groups_to_exclude = project_params[["study_groups_to_exclude"]]
      ),
    packages   =
      c(
        "readr",
        "dplyr",
        "janitor",
        "purrr",
        "forcats",
        "lubridate"
      )
  ),

  tar_target(
    name       = seq_file_directory,
    command    = project_params[["sequencing_file_directory"]],
    format     = "file"
  ),

  tar_target(
    name       = tx_files,
    command    =
      import_counts(
        directory = seq_file_directory,
        metadata  = md
      ),
    packages   =
      c(
        "purrr",
        "magrittr",
        "stringr"
      )
  ),

  tar_target(
    name       = annotation_file,
    command    = project_params[["annotation_file"]],
    format     = "file"
  ),

  tar_target(
    name       = annot,
    command    = read_csv(annotation_file),
    packages   = "readr"
  ),

  tar_target(
    name       = final_md,
    command    =
      create_final_md(
        md               = md,
        tx_files         = tx_files,
        comparison_group = project_params[["comparison_grouping_variable"]],
        control_group    = project_params[["control_group"]]
      ),
    packages   = c(
      "forcats",
      "stringr",
      "dplyr",
      "rlang",
      "tibble"
    )
  ),

  tar_target(
    name       = tx_counts,
    command    =
      import_counts(
        directory = seq_file_directory
      ),
    packages =
      c(
        "tximport",
        "data.table",
        "magrittr",
        "janitor",
        "purrr",
        "stringr"
      )
  ),

  tar_target(
    name = imported_data,
    command = process_counts(
      count_files                  = tx_counts,
      sample_metadata              = final_md,
      study_design                 = project_params[["study_design"]],
      batch_variable               = project_params[["batch_variable"]],
      comparison_grouping_variable = project_params[["comparison_groups"]],
      reference_gene_annotation    = annot,
      aligner                      = project_params[["aligner"]],
      minimum_gene_count           = project_params[["minimum_gene_count"]],
      pc1_zscore_threshold         = project_params[["pc1_zscore_threshold"]],
      pc2_zscore_threshold         = project_params[["pc2_zscore_threshold"]],
      BPPARAM                      = BPPARAM,
      num_sva                      = project_params[["sva_num"]],
      use_combat                   = project_params[["use_combat"]],
      method                       = project_params[["process_method"]]
    ),
    packages =
      c(
        "purrr",
        "tibble",
        "HGNChelper",
        "dplyr",
        "matrixStats",
        "tidyselect",
        "stringr",
        "rlang"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = qc_pca,
    command = imported_data[["qc_pca"]]
  ),

  tar_target(
    name = outlier_samples,
    command = imported_data[["outlier_samples"]]
  ),

  tar_target(
    name       = sva_graph_data,
    command    = plot_sva(imported_data[["sva_graph_data"]])
  ),

  tar_target(
    name = vsc_exprs,
    command =
      imported_data[["variance_stabilized_counts"]] %>%
      as_tibble(rownames = "gene") %>%
      mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
      filter(!is.na(hugo)) %>%
      group_by(hugo) %>%
      slice(1) %>%
      ungroup() %>%
      select(
        -gene,
        gene = hugo
        ) %>%
      column_to_rownames("gene") %>%
      as.matrix()
  ),

  # This should be changed into a list that we can walk through
  tar_target(
    name = banchereau_module_file,
    command = project_params[["banchereau_modules"]],
    format = "file"
  ),

  tar_target(
    name = banchereau_module_annotations_file,
    command = project_params[["banchereau_module_annotations"]],
    format = "file"
  ),

  tar_target(
    name = ldg_module_file,
    command = project_params[["ldg_modules"]],
    format = "file"
  ),

  tar_target(
    name = metasignature_module_file,
    command = project_params[["metasignature_modules"]],
    format = "file"
  ),

  tar_target(
    name = banchereau_modules,
    command = create_module_list(banchereau_module_file)
  ),

  tar_target(
    name = module_annotation,
    command =
      read_csv(banchereau_module_annotations_file) %>%
      mutate(type = as_factor(type))
  ),

  tar_target(
    name = ldg_modules,
    command = create_module_list(ldg_module_file)
  ),

  tar_target(
    name = metasignature_module,
    command = create_module_list(metasignature_module_file)
  ),

  # TODO: can scoreEigengenes handle anything other than a DESeqDataSet or a plain matrix?
  tar_target(
    name = dataset_with_scores,
    command =
      scoreEigengenes(
        object = imported_data[["dataset"]],
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

  # TODO: can extract_module_scores handle anything other than a DESeqDataSet?
  tar_target(
    name = module_scores,
    command =
      extract_module_scores(
        dds = dataset_with_scores,
        names(banchereau_modules),
        names(ldg_modules),
        names(metasignature_module)
      )
  ),

  tar_target(
    name = disp_plot,
    command = plot_dispersion_estimate(dataset_with_scores)
  ),

  tar_target(
    name = annotated_modules,
    command =
      filter(
        .data = module_annotation,
        type != "Undetermined"
      ) %>%
      mutate(type = fct_drop(type))
  ),

  tar_target(
    name = annotated_mod_list,
    command =
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
    name = annotated_module_scores,
    command =
      select(
        .data = module_scores,
        sample_name,
        one_of(annotated_modules$module)
      )
  ),

  tar_target(
    name = module_tbl,
    command =
      create_module_table(
        banchereau_modules,
        ldg_modules,
        metasignature_module,
        module_annotation = module_annotation
      )
  ),

  tar_target(
    name = module_ISGs,
    command =
      extract_module_genes(
        module_table      = module_tbl,
        exprs_mat         = vsc_exprs,
        module_annotation = "Interferon"
      )
  ),

  tar_target(
    name = sample_dists,
    command = parallelDist(t(vsc_exprs)),
    packages = "parallelDist"
  ),

  tar_target(
    name = sample_dendrogram,
    command = as.dendrogram(hclust(sample_dists))
  ),

  tar_target(
    name = sample_cluster_info,
    command =
      ident_clusters(
        irlba(
          A = vsc_exprs,
          nv = 100
          ) %>%
          pluck("v") %>%
          as.data.frame() %>%
          set_colnames(
            paste0(
              "PC",
              seq(100)
              )
            ) %>%
          set_rownames(
            colnames(vsc_exprs)
            ),
        K.max = 20
      )
  ),

  tar_target(
    name = clusters,
    command =
      mutate(
        .data = sample_cluster_info$clusters,
        cluster = as_factor(cluster)
      )
  ),

  # TODO: again, needs to work with things other than a DESeqDataSet
  # TODO: maybe a generalized get metadata from object function?
  tar_target(
    name = study_md,
    command =
      left_join(
        as_tibble(
          colData(dataset_with_scores),
          rownames = "sample_name"),
        clusters
      )
  ),

  tar_target(
    name = annotation_info,
    command =
      column_to_rownames(
        select(
          .data = study_md,
          study_group,
          sex,
          cluster,
          sample_name
        ),
        var = "sample_name"
      )
  ),

  tar_target(
    pca_results,
    run_pca(
      expr_data = vsc_exprs,
      metadata =
        as_tibble(
          x = colData(dataset_with_scores),
          rownames="sample_name"
        ),
      cluster_info = clusters
    )
  ),

  tar_target(
    umap_results,
    run_umap(
      expr_data = vsc_exprs,
      metadata =
        as_tibble(
          x = colData(dataset_with_scores),
          rownames = "sample_name"
        ),
      cluster_info = clusters
    )
  ),

  # tar_target(
  #   name = comparison_results_list,
  #   keep(
  #     resultsNames(dataset_with_scores),
  #     str_detect(
  #       string = resultsNames(dataset_with_scores),
  #       pattern = project_params[["comparison_grouping_variable"]]
  #     )
  #   )
  # ),

  tar_target(
    name = res,
    command =
      create_results_list(
        comparison_list = imported_data[["comparisons"]],
        dds = dataset_with_scores,
        comparison_grouping_variable = project_params[["comparison_groupings"]]
      ),
    cue = tar_cue(mode = "never")
  ),

  # TODO: generalize create_deg_tables?
  tar_target(
    name = down_tables,
    command = create_deg_tables(
      deg_res = imported_data[["res"]],
      comparison_list = imported_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]],
      direction = "down"
    )
  ),

  tar_target(
    name = up_tables,
    command = create_deg_tables(
      deg_res = imported_data[["res"]],
      comparison_list = imported_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]],
      direction = "up"
    )
  ),

  tar_target(
    name = degs,
    command = extract_de_genes(
      results = imported_data[["res"]],
      comparison_list = imported_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]]
    )
  ),

  tar_target(
    name = deg_class,
    command =
      group_degs(
        degs = degs,
        comparison_vars = project_params[["comparison_groups"]]
        )
  ),

  # tar_target(
  #   name = deg_means,
  #   command = calc_deg_means(
  #     exprs = vsc_exprs,
  #     deg_class = deg_class,
  #     metadata = as_tibble(
  #       colData(dataset_with_scores),
  #       rownames = "name"
  #     ) %>%
  #       select(
  #         name,
  #         sex,
  #         ethnicity,
  #         one_of(project_params[["comparison_groups"]])
  #       ),
  #     grouping_variable = study_group
  #   )
  # ),

  tar_target(
    name    = top_degs,
    command =
      extract_top_degs(
        up_tables   = up_tables,
        down_tables = down_tables
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
      rownames(vsc_exprs)
    )
  ),

  tar_target(
    name = vsc_top,
    command = top_variable_genes(
      exprs = vsc_exprs,
      n = 20000
    )
  ),

  tar_target(
    name = sft,
    command = pickSoftThreshold(
      data = vsc_top,
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
  ),

  tar_target(
    name    = wgcna_modules,
    command = blockwiseModules(
      datExpr           = vsc_top,
      power             = find_softPower(sft),
      maxBlockSize      = 20000,
      mergeCutHeight    = 0.2,
      minModuleSize     = 20,
      pamRespectsDendro = FALSE,
      saveTOMs          = FALSE,
      verbose           = 3,
      detectCutHeight   = 0.995,
      TOMDenom          = "min",
      networkType       = "signed hybrid",
      reassignThreshold = 1e-6
    )
  ),

  tar_target(
    name = wgcna_module_genes,
    command =
      enframe(
        x = wgcna_modules$colors,
        name = "gene",
        value = "module"
      )
  ),

  tar_target(
    name = wgcna_module_colors,
    command =
      unique(
        pull(
          .data = wgcna_module_genes,
          module
        )
      )
  ),

  tar_target(
    name = wgcna_hub_genes,
    command =
      chooseTopHubInEachModule(
        datExpr = vsc_top,
        colorh  = wgcna_modules$colors,
        power   = 4,
        type    = "signed hybrid"
      )
  ),

  tar_target(
    name = wgcna_scores,
    command =
      left_join(
        x = select(
          .data =
            as_tibble(
              x        = wgcna_modules$MEs,
              rownames = "sample_name"
            ),
          -MEgrey
        ),
        y =
          as_tibble(
            x        = annotation_info,
            rownames = "sample_name"
          )
      )
  ),

  tar_target(
    name = wgcna_cluster_split,
    command =
      initial_split(
        data = wgcna_scores %>%
          mutate(study_group = fct_drop(study_group)) %>%
          select(
            cluster,
            starts_with("ME")
          ),
        prop = 0.75,
        strata = "cluster"
      )
  ),

  tar_target(
    name = wgcna_cluster_rf_cv,
    command =
      train(
        cluster ~ .,
        method = "parRF",
        data = training(wgcna_cluster_split),
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
  ),

  tar_target(
    name = wgcna_cluster_rf_cv_varImp,
    command =
      varImp(
        object     = wgcna_cluster_rf_cv,
        scale      = FALSE,
        importance = TRUE
      )
  ),

  tar_target(
    name = wgcna_study_group_split,
    command =
      initial_split(
        data =
          select(
            .data         = mutate(
              .data       = wgcna_scores,
              study_group = fct_drop(study_group)
            ),
            study_group,
            starts_with("ME")
          ),
        prop = 0.75,
        strata = "study_group"
      )
  ),

  tar_target(
    name    = wgcna_study_group_train,
    command = training(wgcna_study_group_split)
  ),

  tar_target(
    name    = wgcna_study_group_test,
    command = testing(wgcna_study_group_split)
  ),

  tar_target(
    name    = wgcna_study_group_rf_cv,
    command =
      train(
        form   = study_group ~ .,
        method = "parRF",
        data   = wgcna_study_group_train,
        trControl =
          trainControl(
            method        = "repeatedcv",
            number        = 20,
            repeats       = 20,
            search        = "grid",
            allowParallel = TRUE
          )
      )
  ),

  tar_target(
    name    = wgcna_study_group_rf_cv_varImp,
    command =
      varImp(
        object     = wgcna_study_group_rf_cv,
        scale      = FALSE,
        importance = TRUE
      )
  ),

  tar_target(
    name = filtered_wgcna_module_genes,
    command =
      filter(
        .data =
          mutate(
            .data = wgcna_module_genes,
            hugo  = checkGeneSymbols(gene)[["Suggested.Symbol"]]
          ),
        !is.na(hugo)
      )
  ),

  tar_target(
    name    = MEenriched_list,
    command = module_gsea(
      module_genes       = filtered_wgcna_module_genes,
      module_of_interest = wgcna_module_colors
    ),
    pattern   = map(wgcna_module_colors),
    iteration = "list"
  ),

  tar_target(
    name    = MEenriched,
    command = bind_rows(MEenriched_list)
  ),

  tar_target(
    name    = MEplotting,
    command = module_gsea_plots(enriched_genes = MEenriched)
  ),

  tar_target(
    name    = module_scores_with_md,
    command =
      left_join(
        x = module_scores,
        y = as_tibble(
          x        = annotation_info,
          rownames = "sample_name"
        )
      )
  ),

  tar_target(
    name    = module_cluster_split,
    command =
      initial_split(
        data =
          select(
            .data =
              mutate(
                .data       = module_scores_with_md,
                study_group = fct_drop(study_group)
              ),
            cluster,
            one_of(names(banchereau_modules))
          ),
        prop   = 0.75,
        strata = "cluster"
      )
  ),

  tar_target(
    name    = module_cluster_train,
    command = training(module_cluster_split)
  ),

  tar_target(
    name    = module_cluster_test,
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
            method        = "repeatedcv",
            number        = 10,
            repeats       = 10,
            search        = "grid",
            allowParallel = TRUE
          ),
        importance = TRUE
      )
  ),

  tar_target(
    name = module_cluster_rf_cv_varImp,
    command =
      varImp(
        object     = module_cluster_rf_cv,
        scale      = FALSE,
        importance = TRUE
      )
  ),

  tar_target(
    name = module_study_group_split,
    command =
      initial_split(
        data =
          select(
            .data =
              mutate(
                .data       = module_scores_with_md,
                study_group = fct_drop(study_group)
              ),
            study_group,
            one_of(names(banchereau_modules))
          ),
        prop = 0.75,
        strata = "study_group"
      )
  ),

  tar_target(
    name    = module_study_group_train,
    command = training(module_study_group_split)
  ),

  tar_target(
    name    = module_study_group_test,
    command = testing(module_study_group_split)
  ),

  tar_target(
    name = module_study_group_rf_cv,
    command =
      train(
        form = study_group ~ .,
        method = "parRF",
        data = module_study_group_train,
        trControl =
          trainControl(
            method        = "repeatedcv",
            number        = 10,
            repeats       = 10,
            search        = "grid",
            allowParallel = TRUE
          ),
        importance = TRUE
      )
  ),

  tar_target(
    name    = module_study_group_rf_cv_varImp,
    command =
      varImp(
        object     = module_study_group_rf_cv,
        scale      = FALSE,
        importance = TRUE
      )
  ),

  # tar_target(
  #   name = viral_exprs,
  #   command =
  #     extract_viral_expression(
  #       annotations = annot,
  #       exprs = vsc_exprs,
  #       dds = dataset_with_scores
  #     )
  # ),

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
  ),

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
    name = module_scores_with_viral,
    command =
      reduce(
        .x =
          list(
            module_scores,
            wgcna_scores,
            # viral_exprs,
            clusters
          ),
        .f = left_join
        )
  ),

  tar_target(
    name = annotated_module_scores_with_cluster_class,
    command =
      select(
        .data =
          mutate(
            .data = module_scores_with_viral,
            cluster = as_factor(cluster)
          ),
        cluster,
        study_group,
        one_of(annotated_modules$module)
      )
  ),

  tar_target(
    name = renamed_annotated_module_scores,
    command =
      set_names(
        nm =
          drake_recode(
            target_list             = names(annotated_module_scores_with_cluster_class),
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
        compare_by         = "cluster"
      )
  ),

  tar_target(
    name    = annotated_module_stats_by_disease,
    command =
      modules_compare_with_stats(
        module_score_table = annotated_module_scores_pivot,
        compare_by         = "study_group"
      )
  ),

  tar_target(
    name    = module_scores_pivot,
    command = pivot_module_scores(module_scores = module_scores_with_viral)
  ),

  tar_target(
    name = module_stats_by_cluster,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        compare_by = "cluster"
      )
  ),

  tar_target(
    name    = module_stats_by_disease,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        compare_by         = "study_group"
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
        compare_by         = "cluster"
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
                study_group,
                matches("^ME")
              ),
            -study_group,
            names_to     = "module",
            values_to    = "score"
          ),
        study_group = as_factor(study_group)
      )
  ),

  tar_target(
    name = module_scores_with_viral_by_disease_stats,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_with_viral_by_disease,
        compare_by         = "study_group"
      )
  ),

  tar_target(
    name = group_pal,
    command =
      create_palettes(
        annotated_modules = annotated_modules,
        clusters          = clusters,
        annotation_info   = annotation_info,
        deg_class         = deg_class
      )
  ),

  # TODO: make an optional accessory gene list report
  # tar_target(
  #   name = plafker_gene_list_file,
  #   "references/plafker_gene_list.csv",
  #   format = "file"
  # ),

  # tar_target(
  #   name = plafker_gene_list,
  #   command = create_module_list(plafker_gene_list_file)
  # ),

  # tar_target(
  #   name = pathway_exprs,
  #   command =
  #     map(plafker_gene_list, extract_pathway_exprs, exprs_mat = vsc_exprs, metadata = study_md),
  # ),

  # tar_target(
  #   name = pathway_expr_stats,
  #   command = map(pathway_exprs, calc_gex_stats)
  # ),

  # tar_target(
  #   name = dds_with_pathway_scores,
  #   command =
  #     scoreEigengenes(
  #       object = dataset_with_scores,
  #       module_list = plafker_gene_list,
  #       score_func = 'rsvd'
  #     )
  # ),

  # tar_target(
  #   name = pathway_eigenvalues,
  #   command = extract_module_scores(dds_with_pathway_scores, names(plafker_gene_list))
  # ),

  # tar_target(
  #   name = pathway_eigenvalues_long,
  #   command =
  #     pathway_eigenvalues %>%
  #     pivot_longer(
  #       -sample_name,
  #       names_to      = "pathway",
  #       values_to     = "score"
  #     ) %>%
  #     left_join(
  #       select(
  #         .data       = tar_read(study_md),
  #         sample_name,
  #         study_group
  #       )
  #     )
  #   ),

  # tar_target(
  #   name = pathway_eigenvalues_stats,
  #   command =
  #     pathway_eigenvalues_long %>%
  #     group_by(pathway) %>%
  #     wilcox_test(
  #       score ~ study_group,
  #       ref.group = "control"
  #     ) %>%
  #     adjust_pvalue(method = "BH") %>%
  #     add_significance() %>%
  #     grouped_add_xy_positions(
  #       stats_tbl = .,
  #       data_tbl = pathway_eigenvalues_long,
  #       group_var = pathway,
  #       compare_value = score,
  #       percent_shift_down = 0.99
  #     )
  # ),

  # This... did not work out so great.  ICA not the best for figuring out the modules
  # but this could be adapted for PLIER or NMF or whatever
  # tar_target(
  #   name = ica_res,
  #   command =
  #     icafast(
  #       X = vsc_top,
  #       nc = 25,
  #       alg = "par"
  #       ),
  #   packages = "ica"
  # ),
  #
  # tar_target(
  #   name = ica_eigenvalues,
  #   command =
  #     pluck(
  #       .x = ica_res,
  #       "S"
  #     ) %>%
  #     as.data.frame() %>%
  #     set_rownames(rownames(vsc_top)) %>%
  #     set_colnames(str_glue("IC{1:25}"))
  # ),
  #
  # tar_target(
  #   name = ica_eigenvalues_long,
  #   command =
  #     as_tibble(
  #       x = ica_eigenvalues,
  #       rownames = "sample_name") %>%
  #     pivot_longer(
  #       cols = starts_with("IC"),
  #       names_to = "IC",
  #       values_to = "score"
  #       ) %>%
  #     left_join(
  #       y =
  #         select(
  #           .data = study_md,
  #           sample_name,
  #           cluster,
  #           study_group
  #           )
  #       )
  # ),
  #
  # tar_target(
  #   name = ica_eigenvalues_stats,
  #   command =
  #     group_by(
  #       .data = ica_eigenvalues_long,
  #       pathway
  #       ) %>%
  #     wilcox_test(
  #       score ~ study_group,
  #       ref.group = "control"
  #       ) %>%
  #     adjust_pvalue(method = "BH") %>%
  #     add_significance() %>%
  #     grouped_add_xy_positions(
  #       stats_tbl = .,
  #       data_tbl = ica_eigenvalues_long,
  #       group_var = IC,
  #       compare_value = score,
  #       percent_shift_down = 0.99
  #       )
  #   ),
  #
  # tar_target(
  #   name = ica_loadings,
  #   command =
  #     pluck(
  #       .x = ica_res,
  #       "M"
  #       ) %>%
  #     as.data.frame() %>%
  #     set_rownames(colnames(vsc_top)) %>%
  #     set_colnames(str_glue("IC{1:25}"))
  #   ),
  #
  # tar_target(
  #   name = ic_top_25_genes,
  #   command =
  #     ica_loadings %>%
  #     as_tibble(rownames = "gene") %>%
  #     pivot_longer(
  #       -gene,
  #       names_to = "IC",
  #       values_to = "score"
  #     ) %>%
  #     group_by(IC) %>%
  #     top_n(
  #       n = 25,
  #       wt = score
  #     ) %>%
  #     select(-score) %>%
  #     rename(module = "IC") %>%
  #     mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
  #     filter(!is.na(hugo))
  #   ),
  #
  # tar_target(
  #   name = ic_enriched_list,
  #   command =
  #     map_dfr(
  #       .x = unique(ic_top_25_genes$module),
  #       .f = function(x){
  #         module_gsea(
  #           module_genes = ic_top_25_genes,
  #           module_of_interest = x)
  #       }
  #     ),
  #   packages = "purrr"
  #   ),
  #
  # tar_target(
  #   name = ic_plotting,
  #   command =
  #     module_gsea_plots(ic_enriched_list) %>%
  #     mutate(
  #       module =
  #         str_remove(
  #           string = module,
  #           pattern = "^ME"
  #           )
  #       )
  # ),

  tar_target(
    name     = output_expression,
    command  =
      save_table_to_disk(
        file_to_output =
          as_tibble(
            x = vsc_exprs,
            rownames = "sample_name"
          ),
        output_name    = "processed_data/variance_stabilized_expression.csv.gz"
        ),
    format   = "file"
  ),

  tar_target(
    name     = output_metadata,
    command  =
      save_table_to_disk(
        file_to_output =
          as_tibble(
            x        = colData(dataset_with_scores),
            rownames = "sample_name"
          ),
        output_name    = "processed_data/sample_metadata.csv.gz"
      ),
    format   = "file"
  ),

  tar_target(
    name     = output_module_scores,
    command  =
      save_table_to_disk(
        file_to_output = module_scores,
        output_name    = "processed_data/module_scores.csv.gz"
      ),
    format   = "file"
  ),

  tar_target(
    name     = output_wgcna_scores,
    command  =
      save_table_to_disk(
        file_to_output = wgcna_scores,
        output_name    = "processed_data/wgcna_scores.csv.gz"
      ),
    format   = "file"
  ),

  # tar_target(
  #   name     = output_pathway_eigenvalues,
  #   command  =
  #     save_table_to_disk(
  #       file_to_output = pathway_eigenvalues,
  #       output_name    = "processed_data/plafker_pathway_scores.csv.gz"
  #     ),
  #   format   = "file"
  # ),

  tar_target(
    name     = output_wgcna_module_genes,
    command  =
      save_table_to_disk(
        file_to_output = wgcna_module_genes,
        output_name    = "processed_data/wgcna_module_genes.csv.gz"
      ),
    format   = "file"
  ),

  tar_render(
    name          = primary_report,
    path          = "analysis/report.rmd",
    params        = list(
      set_title   = "Initial COVID PCV samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
    output_dir    = "reports/"
  ),

  tar_render(
    name          = qc_report,
    path          =  "analysis/qc_report.rmd",
    params        = list(
      set_title   = "Initial COVID PCV samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
    output_dir    = "reports/"
  )

  # tar_render(
  #   name          = pathway_report,
  #   path          =  "analysis/pathway_gex_plots.rmd",
  #   params        = list(
  #     set_title   = "Expression of pathway genes",
  #     set_author  = "Miles Smith"
  #   ),
  #   output_dir    = "reports/",
  #   packages      = c("tidyverse", "markdown", "knitr", "ComplexHeatmap", "targets", "grid", "magrittr", "rlang", "viridis", "circlize", "ggpubr", "ggbeeswarm", "paletteer", "ggtext", "tidyHeatmap")
  # )
)
