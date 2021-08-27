library(targets)
library(tarchetypes)

source("code/project_parameters.R")

source("code/plan/01_import_funcs.R")
source("code/plan/02_filtering_funcs.R")
source("code/plan/04_module_funcs.R")
source("code/plan/06_dimensional_reduction_funcs.R")
source("code/plan/07_differential_expression_funcs.R")
source("code/plan/08_WGCNA_funcs.R")
source("code/plan/09_ml_funcs.R")
source("code/plan/10_viral_transcript_funcs.R")
source("code/plan/11_stats_testing_funcs.R")
source("code/plan/13_pathways.R")
source("code/plan/97_misc_functions.R")
source("code/plan/98_palettes_funcs.R")
source("code/plan/99_output_funcs.R")

options(tidyverse.quiet = TRUE)
future::plan(strategy = future::multisession)

list(
  tar_target(
    name = raw_metadata,
    command    = project_params[["metadata_file"]],
    format     = "file",
    deployment = "main",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = md,
    command =
      import_metadata(
        metadata_file                = raw_metadata,
        comparison_grouping_variable = project_params[["comparison_grouping_variable"]],
        projects_to_include          = project_params[["projects_to_include"]],
        projects_to_exclude          = project_params[["projects_to_exclude"]],
        groups_to_include            = project_params[["groups_to_include"]],
        groups_to_exclude            = project_params[["groups_to_exclude"]],
        sample_name_column           = project_params[["sample_name_column"]],
        grouping_column              = project_params[["grouping_column"]],
        project_column               = project_params[["project_column"]],
        regression_columns           = project_params[["regression_columns"]],
        filter_column                = project_params[["filter_column"]],
        samples_to_exclude           = project_params[["manual_sample_removal"]],
        extra_controls_metadata_file = raw_metadata
      ),
    packages =
      c(
        "readr",
        "readxl",
        "dplyr",
        "janitor",
        "purrr",
        "forcats",
        "stringr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = seq_file_directory,
    command    = project_params[["sequencing_file_directory"]],
    format     = "file",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = tx_files,
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
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = annotation_file,
    command    = project_params[["annotation_file"]],
    format     = "file",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = annot,
    command    = read_csv(annotation_file),
    packages   = "readr",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = final_md,
    command    =
      create_final_md(
        md               = md,
        tx_files         = tx_files,
        comparison_group = project_params[["comparison_grouping_variable"]],
        control_group    = project_params[["control_group"]],
        sample_name = project_params[["sample_name_column"]]
      ),
    packages = c(
      "forcats",
      "stringr",
      "dplyr",
      "rlang",
      "tibble"
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = imported_data,
    command =
      prep_data_import(
        count_files        = tx_files,
        sample_metadata    = final_md,
        aligner            = project_params[["aligner"]],
        annotations        = annot,
        minimum_gene_count = 1,
        removal_pattern    = "^RNA5",
        only_hugo          = project_params[["only_hugo_named_genes"]]
      ),
    cue = tar_cue(mode = "never"),
    packages =
      c(
        "magrittr",
        "tximport",
        "data.table",
        "glue",
        "tibble",
        "dplyr",
        "purrr",
        "HGNChelper",
        "stringr"
      )
  ),

  tar_target(
    name = processed_data,
    command = process_counts(
      imported_counts              = imported_data,
      batch_variable               = project_params[["batch_variable"]],
      comparison_grouping_variable = project_params[["grouping_column"]],
      study_design                 = project_params[["study_design"]],
      pc1_zscore_threshold         = project_params[["pc1_zscore_threshold"]],
      pc2_zscore_threshold         = project_params[["pc2_zscore_threshold"]],
      BPPARAM                      = BPPARAM,
      use_combat                   = project_params[["use_combat"]],
      minimum_gene_count           = project_params[["minimum_gene_count"]],
      num_sva                      = project_params[["sva_num"]],
      method                       = project_params[["process_method"]],
      control_group                = project_params[["control_group"]]
    ),
    packages =
      c(
        "edgeR",
        "DESeq2",
        "sva",
        "Rfast",
        "limma",
        "dplyr",
        "purrr",
        "tibble",
        "stringr",
        "rlang",
        "tidyr",
        "magrittr",
        "gtools"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = qc_pca,
    command = processed_data[["qc_pca"]],
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = outlier_samples,
    command = processed_data[["outlier_samples"]],
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = sva_graph_data,
    command = plot_sva(processed_data[["sva_graph_data"]]),
    packages =
      c(
        "ggplot2",
        "cowplot",
        "tidyselect",
        "tidyr",
        "tibble",
        "purrr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = vsc_exprs,
    command =
      extract_transformed_data(
        data_obj = processed_data[["variance_stabilized_counts"]]
        ),
    packages =
      c(
        "tibble",
        "dplyr",
        "HGNChelper"
      ),
    cue = tar_cue(mode = "never")
  ),

  # This should be changed into a list that we can walk through
  tar_target(
    name = banchereau_module_file,
    command = project_params[["banchereau_modules"]],
    format = "file",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = banchereau_module_annotations_file,
    command = project_params[["banchereau_module_annotations"]],
    format = "file",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = ldg_module_file,
    command = project_params[["ldg_modules"]],
    format = "file",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = metasignature_module_file,
    command = project_params[["metasignature_modules"]],
    format = "file",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = banchereau_modules,
    command = create_module_list(banchereau_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = module_annotation,
    command =
      readr::read_csv(banchereau_module_annotations_file) %>%
      dplyr::mutate(type = forcats::as_factor(type)),
    packages =
      c(
        "readr",
        "dplyr",
        "forcats",
        "magrittr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = ldg_modules,
    command = create_module_list(ldg_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = metasignature_module,
    command = create_module_list(metasignature_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      ),
    cue = tar_cue(mode = "never")
  ),


  tar_target(
    name = dataset_with_scores,
    command =
      scoreEigengenes(
        object = processed_data[["dataset"]],
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
      ),
    packages =
      c(
        "moduleScoreR",
        "magrittr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = module_scores,
    command =
      extract_module_scores(
        object = dataset_with_scores,
        names(banchereau_modules),
        names(ldg_modules),
        names(metasignature_module)
      ) %>%
      left_join(
        as_tibble(
          annotation_info,
          rownames = "sample_name"
          )
      ),
    packages =
      c(
        "purrr",
        "rlang",
        "dplyr",
        "tibble",
        "tidyselect",
        "magrittr"
      ),
    cue = tar_cue(mode = "never")
  ),

  # tar_target(
  #   name = disp_plot,
  #   command = plot_dispersion_estimate(dataset_with_scores),
  #   cue = tar_cue(mode = "never")
  # ),

  tar_target(
    name = annotated_modules,
    command =
      dplyr::filter(
        .data = module_annotation,
        type != "Undetermined"
      ) %>%
      dplyr::mutate(type = forcats::fct_drop(type)),
    packages =
      c(
        "dplyr",
        "forcats",
        "magrittr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = annotated_mod_list,
    command = {
      dplyr::mutate(
        .data = annotated_modules,
        module_type =
          paste(
            module,
            type,
            sep = " - "
          )
      ) %>%
      dplyr::select(-type) %>%
      tibble::deframe()
      },
    packages =
      c(
        "dplyr",
        "tibble",
        "magrittr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = annotated_module_scores,
    command =
      dplyr::select(
        .data = module_scores,
        sample_name,
        tidyselect::one_of(annotated_modules$module)
      ),
    packages =
      c(
        "dplyr",
        "tidyselect"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = module_tbl,
    command =
      create_module_table(
        banchereau_modules,
        ldg_modules,
        metasignature_module,
        module_annotation = module_annotation
      ),
    packages =
      c(
        "rlang",
        "tibble",
        "purrr",
        "tidyr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = module_ISGs,
    command =
      extract_module_genes(
        module_table      = module_tbl,
        exprs_mat         = vsc_exprs,
        module_annotation = "Interferon"
      ),
    packages =
      c(
        "dplyr",
        "tibble"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = sample_dists,
    command = parallelDist(t(vsc_exprs)),
    packages = "parallelDist",
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = sample_dendrogram,
    command = as.dendrogram(hclust(sample_dists)),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = sample_cluster_info,
    command =
      ident_clusters(
        irlba::irlba(
          A = vsc_exprs,
          nv = 100
        ) %>%
          purrr::chuck("v") %>%
          as.data.frame() %>%
          magrittr::set_colnames(
            paste0(
              "PC",
              seq(100)
            )
          ) %>%
          magrittr::set_rownames(
            colnames(vsc_exprs)
          ),
        K.max = 20
      ),
    packages =
      c(
        "randomForest",
        "cluster",
        "stats",
        "tibble",
        "forcats",
        "dplyr",
        "magrittr"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = study_md,
    command =
      dplyr::left_join(
        getMetaData(dataset_with_scores),
        sample_cluster_info[["clusters"]]
      ),
    packages = c(
      "dplyr",
      "SummarizedExperiment",
      "tibble"
    ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = annotation_info,
    command = {
      tibble::column_to_rownames(
        dplyr::select(
          .data = study_md,
          tidyselect::one_of(
            project_params[["grouping_column"]],
            project_params[["project_column"]],
            project_params[["regression_columns"]]
          ),
          sample_name
        ),
        var = "sample_name"
      )
    },
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tibble"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    pca_results,
    run_pca(
      expr_data    = vsc_exprs,
      metadata     = study_md,
      cluster_info = clusters
    ),
    packages =
      c(
        "irlba",
        "tibble",
        "purrr",
        "dplyr"
      )
  ),

  tar_target(
    umap_results,
    run_umap(
      expr_data    = vsc_exprs,
      metadata     = study_md,
      cluster_info = clusters
    ),
    packages =
      c(
        "uwot",
        "tibble",
        "parallel",
        "magrittr",
        "dplyr"
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

  # TODO: {targets} should be able to handle mapping
  # the comparison list to the function instead of us mapping it
  tar_target(
    name = res,
    command =
      create_results_list(
        comparison_list              = processed_data[["comparisons"]],
        object                       = dataset_with_scores,
        comparison_grouping_variable = project_params[["comparison_groupings"]]
      ),
    packages =
      c(
        "purrr",
        "DESeq2",
        "tibble",
        "rlang",
        "magrittr",
        "edgeR",
        "matrixStats",
        "rstatix",
        "dplyr"
      ),
    cue = tar_cue(mode = "never")
  ),

  # TODO: generalize create_deg_tables?
  tar_target(
    name = down_tables,
    command = create_deg_tables(
      deg_res           = processed_data[["res"]],
      comparison_list   = processed_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]],
      direction         = "down"
    ),
    packages =
      c(
        "purrr",
        "tibble",
        "dplyr",
        "tidyselect",
        "rlang",
        "stringr"
      )
  ),

  tar_target(
    name = up_tables,
    command  = create_deg_tables(
      deg_res           = processed_data[["res"]],
      comparison_list   = processed_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]],
      direction         = "up"
    ),
    packages =
      c(
        "purrr",
        "tibble",
        "dplyr",
        "tidyselect",
        "rlang",
        "stringr"
      )
  ),

  tar_target(
    name = degs,
    command  = extract_de_genes(
      results           = processed_data[["res"]],
      comparison_list   = processed_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]]
    ),
    packages =
      c(
        "purrr",
        "tibble",
        "filter",
        "pull",
        "rlang"
      )
  ),

  tar_target(
    name = deg_class,
    command  =
      group_degs(
        degs            = degs,
        comparison_vars = project_params[["comparison_groups"]]
      ),
    packages =
      c(
        "tibble",
        "tidyr",
        "dplyr",
        "stringr"
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
    name = top_degs,
    command  =
      extract_top_degs(
        up_tables   = up_tables,
        down_tables = down_tables
      ),
    packages =
      c(
        "purrr",
        "dplyr",
        "rlang"
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
    command  = top_variable_genes(
      exprs = vsc_exprs,
      n     = 20000
    ),
    packages =
      c(
        "matrixStats",
        "rlang",
        "tibble",
        "dplyr"
      )
  ),

  tar_target(
    name = sft,
    command = WGCNA::pickSoftThreshold(
      data = vsc_top,
      powerVector =
        c(
          seq(10),
          seq(
            from = 12,
            to = 30,
            by = 1
            )
        ),
      verbose = 5
    ),
    packages = "WGCNA"
  ),

  tar_target(
    name = wgcna_modules,
    command = WGCNA::blockwiseModules(
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
    ),
    packages = "WGCNA"
  ),

  tar_target(
    name = wgcna_module_genes,
    command  =
      tibble::enframe(
        x     = wgcna_modules$colors,
        name = "gene",
        value = "module"
      ),
    packages = "tibble"
  ),

  tar_target(
    name = wgcna_module_colors,
    command =
      unique(
        dplyr::pull(
          .data = wgcna_module_genes,
          module
        )
      ),
    packages = "dplyr"
  ),

  tar_target(
    name = wgcna_hub_genes,
    command =
      WGCNA::chooseTopHubInEachModule(
        datExpr = vsc_top,
        colorh  = wgcna_modules$colors,
        power   = 4,
        type    = "signed hybrid"
      ),
    packages = "WGCNA"
  ),

  tar_target(
    name = wgcna_scores,
    command  =
      dplyr::left_join(
        x = dplyr::select(
          .data =
            tibble::as_tibble(
              x        = wgcna_modules$MEs,
              rownames = "sample_name"
            ),
          -MEgrey
        ),
        y =
          tibble::as_tibble(
            x        = annotation_info,
            rownames = "sample_name"
          )
      ),
    packages =
      c(
        "dplyr",
        "tibble"
      )
  ),

  # TODO: can we vectorize the rf model building so that these could
  # be run in parallel?

  tar_target(
    name = cluster_wgcna_rf_model,
    command  =
      rf_classifier(
        dataset            = wgcna_scores,
        classification_var = "cluster",
        starts_with("ME"),
        train_proportion   = 0.75,
        print              = FALSE
      ),
    packages =
      c(
        "rlang",
        "ranger",
        "randomForest",
        "dplyr",
        "forcats",
        "rsample",
        "recipes",
        "parsnip",
        "dials",
        "workflows",
        "tune"
      ),
    cue = tar_cue(mode = "never")
  ),

  tar_target(
    name = comparison_wgcna_rf_model,
    command  =
      rf_classifier(
        dataset            = wgcna_scores,
        classification_var = project_params[["comparison_grouping_variable"]],
        starts_with("ME"),
        train_proportion   = 0.75,
        print              = FALSE
      ),
    packages =
      c(
        "rlang",
        "ranger",
        "randomForest",
        "dplyr",
        "forcats",
        "rsample",
        "recipes",
        "parsnip",
        "dials",
        "workflows",
        "tune"
      ),
    cue      = tar_cue(mode = "never")
  ),


  tar_target(
    name = filtered_wgcna_module_genes,
    command =
      dplyr::filter(
        .data =
          dplyr::mutate(
            .data = wgcna_module_genes,
            hugo  = HGNChelper::checkGeneSymbols(gene)[["Suggested.Symbol"]]
          ),
        !is.na(hugo)
      ),
    packages =
      c(
        "dplyr",
        "HGNChelper"
      )
  ),

  tar_target(
    name = MEenriched_list,
    command = module_gsea(
      module_genes       = filtered_wgcna_module_genes,
      module_of_interest = wgcna_module_colors
    ),
    pattern   = map(wgcna_module_colors),
    iteration = "list",
    packages  =
      c(
        "dplyr",
        "clusterProfiler"
      )
  ),

  tar_target(
    name = MEenriched,
    command  = dplyr::bind_rows(MEenriched_list),
    packages = "dplyr"
  ),

  tar_target(
    name = MEplotting,
    command  = module_gsea_plots(enriched_genes = MEenriched),
    packages =
      c(
        "dplyr",
        "purrr",
        "stringr",
        "magrittr"
      )
  ),

  tar_target(
    name = module_scores_with_md,
    command =
      dplyr::left_join(
        x = module_scores,
        y = tibble::as_tibble(
          x        = annotation_info,
          rownames = "sample_name"
        )
      ),
    packages =
      c(
        "dplyr",
        "tibble"
        )
  ),

  tar_target(
    name = cluster_module_rf_model,
    command  =
      rf_classifier(
        dataset            = module_scores,
        classification_var = "cluster",
        starts_with("ME"),
        train_proportion   = 0.75,
        print              = FALSE
      ),
    packages =
      c(
        "rlang",
        "ranger",
        "randomForest",
        "dplyr",
        "forcats",
        "rsample",
        "recipes",
        "parsnip",
        "dials",
        "workflows",
        "tune"
      ),
    cue      = tar_cue(mode = "never")
  ),

  tar_target(
    name = comparison_module_rf_model,
    command  =
      rf_classifier(
        dataset            = module_scores,
        classification_var = project_params[["comparison_grouping_variable"]],
        starts_with("M"),
        train_proportion   = 0.75,
        print              = FALSE
      ),
    packages =
      c(
        "rlang",
        "ranger",
        "randomForest",
        "dplyr",
        "forcats",
        "rsample",
        "recipes",
        "parsnip",
        "dials",
        "workflows",
        "tune"
      ),
    cue      = tar_cue(mode = "never")
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
      dplyr::pull(
        .data =
          dplyr::filter(
            .data = annotated_modules,
            type == "Interferon"
          ),
        var = module
      ),
    packages = "dplyr"
  ),

  tar_target(
    name = inflame_modules,
    command =
      dplyr::pull(
        .data =
          dplyr::filter(
            .data = annotated_modules,
            type == "Inflammation"
          ),
        var = module
      ),
    packages = "dplyr"
  ),

  tar_target(
    name = module_scores_with_viral,
    command =
      purrr::reduce(
        .x =
          list(
            module_scores,
            wgcna_scores,
            # viral_exprs,
            clusters
          ),
        .f = dplyr::left_join
      ),
    packages =
      c(
        "dplyr",
        "purrr"
      )
  ),

  tar_target(
    name = annotated_module_scores_with_cluster_class,
    command =
      dplyr::select(
        .data = module_scores_with_viral,
        cluster,
        study_group,
        tidyselect::one_of(annotated_modules$module)
      ),
    packages =
      c(
        "dplyr",
        "tidyselect"
      )
  ),

  tar_target(
    name = renamed_annotated_module_scores,
    command =
      rlang::set_names(
        nm =
          targets_recode(
            target_list = names(annotated_module_scores_with_cluster_class),
            thing_to_unquote_splice = annotated_mod_list
            ),
        x  = annotated_module_scores_with_cluster_class
      ),
    packages =
      c(
        "rlang",
        "dplyr"
      )
  ),

  tar_target(
    name = annotated_module_scores_pivot,
    command =
      tidyr::pivot_longer(
        data      = renamed_annotated_module_scores,
        cols      = tidyselect::starts_with("M"),
        names_to  = "module",
        values_to = "score"
      ),
    packages =
      c(
        "tidyr",
        "tidyselect",
        "purrr",
        "rstatix"
      )
  ),

  tar_target(
    name = annotated_module_stats_by_cluster,
    command =
      modules_compare_with_stats(
        module_score_table = annotated_module_scores_pivot,
        compare_by         = "cluster"
      ),
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      )
  ),

  tar_target(
    name = annotated_module_stats_by_comparison_var,
    command =
      modules_compare_with_stats(
        module_score_table = annotated_module_scores_pivot,
        compare_by         = project_params[["comparison_grouping_variable"]]
      ),
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      )
  ),

  tar_target(
    name = module_stats_by_cluster,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        compare_by = "cluster"
      ),
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      )
  ),

  tar_target(
    name = module_stats_by_comparison,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        compare_by         = project_params[["comparison_grouping_variable"]]
      ),
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      )
  ),

  tar_target(
    name = module_scores_with_viral_by_cluster_stats,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_with_viral_by_cluster,
        compare_by         = "cluster"
      ),
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      )
  ),

  tar_target(
    name = module_scores_with_viral_by_comparison_stats,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_with_viral_by_comparison,
        compare_by         = project_params[["comparison_grouping_variable"]]
      )
  ),

  tar_target(
    name = module_scores_pivot,
    command =
      pivot_module_scores(
        module_scores = module_scores_with_viral,
        tidyselect::matches("^M[[:digit:]]+"),
        mg,
        tidy_select::starts_with("ldg")
        ),
    packages =
      c(
        "dplyr",
        "forcats",
        "rlang"
      )
  ),

  tar_target(
    name = wgcna_scores_with_viral_by_cluster,
    command =
      tidyr::pivot_longer(
        data =
          dplyr::select(
            .data = module_scores_with_viral,
            cluster,
            tidyselect::matches("^ME")
          ),
        -cluster,
        names_to     = "module",
        values_to    = "score"
      ),
    packages =
      c(
        "tidyr",
        "tidyselect",
        "dplyr"
      )
  ),

  tar_target(
    name = module_scores_with_viral_by_comparison,
    command =
      tidyr::pivot_longer(
        data =
          dplyr::select(
            .data = module_scores_with_viral,
            project_params[["comparison_grouping_variable"]],
            tidyselect::matches("^ME")
          ),
        -study_group,
        names_to     = "module",
        values_to    = "score"
      ),
    packages =
      c(
        "tidyr",
        "tidyselect",
        "dplyr"
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
    name = output_expression,
    command  =
      writeData(
        object =
          tibble::as_tibble(
            x = vsc_exprs,
            rownames = "sample_name"
          ),
        output_name = "processed_data/variance_stabilized_expression.csv.gz"
      ),
    format   = "file",
    packages =
      c(
        "data.table",
        "tibble"
      )
  ),

  tar_target(
    name = output_metadata,
    command  =
      writeMetaData(
        object = dataset_with_scores,
        output_name = "processed_data/sample_metadata.csv.gz"
      ),
    format   = "file",
    packages =
      c(
        "data.table",
        "pluck"
      ),
  ),

  tar_target(
    name = output_module_scores,
    command  =
      save_table_to_disk(
        file_to_output = module_scores,
        output_name = "processed_data/module_scores.csv.gz"
      ),
    packages = "data.table",
    format   = "file"
  ),

  tar_target(
    name = output_wgcna_scores,
    command  =
      writeData(
        object = wgcna_scores,
        output_name = "processed_data/wgcna_scores.csv.gz"
      ),
    packages = "data.table",
    format   = "file"
  ),

  # tar_target(
  #   name = output_pathway_eigenvalues,
  #   command  =
  #     save_table_to_disk(
  #       file_to_output = pathway_eigenvalues,
  #       output_name = "processed_data/plafker_pathway_scores.csv.gz"
  #     ),
  #   format   = "file"
  # ),

  tar_target(
    name = output_wgcna_module_genes,
    command  =
      writeData(
        object = wgcna_module_genes,
        output_name = "processed_data/wgcna_module_genes.csv.gz"
      ),
    format   = "file",
    packages = "data.table"
  ),

  tar_render(
    name = primary_report,
    path          = "analysis/report.rmd",
    params        = list(
      set_title   = "Initial COVID PCV samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
    output_dir    = "reports/"
  ),

  tar_render(
    name = qc_report,
    path          =  "analysis/qc_report.rmd",
    params        = list(
      set_title   = "Initial COVID PCV samples RNAseq Analysis",
      set_author  = "Miles Smith"
    ),
    output_dir    = "reports/"
  )

  # tar_render(
  #   name = pathway_report,
  #   path          =  "analysis/pathway_gex_plots.rmd",
  #   params        = list(
  #     set_title   = "Expression of pathway genes",
  #     set_author  = "Miles Smith"
  #   ),
  #   output_dir    = "reports/",
  #   packages      = c("tidyverse", "markdown", "knitr", "ComplexHeatmap", "targets", "grid", "magrittr", "rlang", "viridis", "circlize", "ggpubr", "ggbeeswarm", "paletteer", "ggtext", "tidyHeatmap")
  # )
)
