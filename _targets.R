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
# source("code/plan/12_table_outputs.R")
source("code/plan/13_pathways.R")
source("code/plan/14_report_funcs.R")
source("code/plan/97_misc_functions.R")
source("code/plan/98_palettes_funcs.R")
source("code/plan/99_output_funcs.R")

options(tidyverse.quiet = TRUE)
future::plan(strategy = future::multisession)
utils::globalVariables("where")

list(

  targets::tar_target(
    name = comparison_groups,
    command = project_params[["comparison_groups"]]
  ),

  targets::tar_target(
    name = sample_species_org,
    command = findOrgDb(project_params[["sample_species"]]),
  ),

  targets::tar_target(
    name = raw_metadata,
    command    = project_params[["metadata_file"]],
    format     = "file",
    deployment = "main"
  ),

  targets::tar_target(
    name = md,
    command =
      import_metadata(
        metadata_file                 = raw_metadata,
        metadata_sheet                = project_params[["metadata_sheet"]],
        comparison_grouping_variable  = project_params[["comparison_grouping_variable"]],
        projects_to_include           = project_params[["projects_to_include"]],
        projects_to_exclude           = project_params[["projects_to_exclude"]],
        groups_to_include             = project_params[["groups_to_include"]],
        groups_to_exclude             = project_params[["groups_to_exclude"]],
        sample_name_column            = project_params[["sample_name_column"]],
        grouping_column               = project_params[["grouping_column"]],
        project_column                = project_params[["project_column"]],
        regression_columns            = project_params[["regression_columns"]],
        filter_column                 = project_params[["filter_column"]],
        filter_value                  = project_params[["filter_value"]],
        extra_columns                 = project_params[["extra_columns"]],
        samples_to_exclude            = project_params[["manual_sample_removal"]],
        extra_controls_metadata_file  = raw_metadata,
        extra_controls_metadata_sheet = project_params[["main_sample_sheet"]],
        extra_controls_metadata_skip  = project_params[["main_sample_sheet_skip"]],
        skip_lines                    = project_params[["skip_lines"]],
        control_group                 = project_params[["control_group"]]
      ),
    packages =
      c(
        "readr",
        "readxl",
        "dplyr",
        "purrr",
        "forcats",
        "stringr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = seq_file_directory,
    command    = project_params[["sequencing_file_directory"]],
    format     = "file"
  ),

  targets::tar_target(
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
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = annotation_file,
    command    = project_params[["annotation_file"]],
    format     = "file"
  ),

  targets::tar_target(
    name = annot,
    command    = read_csv(annotation_file),
    packages   = "readr"
  ),

  targets::tar_target(
    name = final_md,
    command    =
      create_final_md(
        md               = md,
        tx_files         = tx_files,
        comparison_group = project_params[["comparison_grouping_variable"]],
        control_group    = project_params[["control_group"]],
        sample_name      = project_params[["sample_name_column"]]
      ),
    packages = c(
      "forcats",
      "stringr",
      "dplyr",
      "rlang",
      "tibble"
    ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
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
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
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
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = qc_pca,
    command = processed_data[["qc_pca"]]
  ),

  targets::tar_target(
    name = outlier_samples,
    command = processed_data[["outlier_samples"]]
  ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
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
      )
  ),

  # This should be changed into a list that we can walk through
  targets::tar_target(
    name = banchereau_module_file,
    command = project_params[["banchereau_modules"]],
    format = "file"
  ),

  targets::tar_target(
    name = banchereau_module_annotations_file,
    command = project_params[["banchereau_module_annotations"]],
    format = "file"
  ),

  targets::tar_target(
    name = ldg_module_file,
    command = project_params[["ldg_modules"]],
    format = "file"
  ),

  targets::tar_target(
    name = metasignature_module_file,
    command = project_params[["metasignature_modules"]],
    format = "file"
  ),

  targets::tar_target(
    name = banchereau_modules,
    command = create_module_list(banchereau_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      )
  ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
    name = ldg_modules,
    command = create_module_list(ldg_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      )
  ),

  targets::tar_target(
    name = metasignature_module,
    command = create_module_list(metasignature_module_file),
    packages =
      c(
        "purrr",
        "tibble",
        "tidyr",
        "readr",
        "dplyr"
      )
  ),


  targets::tar_target(
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
      )
  ),

  targets::tar_target(
    name = module_scores,
    command =
      extract_module_scores(
        object = dataset_with_scores,
        names(banchereau_modules),
        names(ldg_modules),
        names(metasignature_module)
      ),
    packages =
      c(
        "purrr",
        "rlang",
        "dplyr",
        "tibble",
        "tidyselect",
        "magrittr"
      )
  ),

  # targets::tar_target(
  #   name = disp_plot,
  #   command = plot_dispersion_estimate(dataset_with_scores),
  #   cue = tar_cue(mode = "never")
  # ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
    name = annotated_module_scores,
    command =
      dplyr::select(
        .data = module_scores,
        sample_name,
        tidyselect::one_of(annotated_modules[["module"]])
      ),
    packages =
      c(
        "dplyr",
        "tidyselect"
      )
  ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
    name = sample_dists,
    command = parallelDist(t(vsc_exprs)),
    packages = "parallelDist"
  ),

  targets::tar_target(
    name = sample_dendrogram,
    command = as.dendrogram(hclust(sample_dists))
  ),

  targets::tar_target(
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
      )
  ),

  targets::tar_target(
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
    )
  ),

  targets::tar_target(
    name = annotation_info,
    command = {
      tibble::column_to_rownames(
        dplyr::select(
          .data = study_md,
          tidyselect::one_of(
            project_params[["grouping_column"]],
            project_params[["project_column"]],
            project_params[["regression_columns"]],
            "cluster"
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
      )
  ),

  targets::tar_target(
    pca_results,
    run_pca(
      expr_data    = vsc_exprs,
      metadata     = study_md
    ),
    packages =
      c(
        "irlba",
        "tibble",
        "purrr",
        "dplyr"
      )
  ),

  targets::tar_target(
    umap_results,
    run_umap(
      expr_data    = vsc_exprs,
      metadata     = study_md
    ),
    packages =
      c(
        "uwot",
        "tibble",
        "parallel",
        "rlang",
        "dplyr"
      )
  ),

  # TODO: {targets} should be able to handle mapping
  # the comparison list to the function instead of us mapping it
  targets::tar_target(
    name = res,
    command =
      create_results_list(
        comparison_list              = processed_data[["comparisons"]],
        method                       = project_params[['process_method']],
        object                       = dataset_with_scores,
        comparison_grouping_variable = project_params[["comparison_grouping_variable"]],
        BPPARAM                      = project_params[["BPPARAM"]],
        design                       = processed_data[['design_matrix']]
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
        "dplyr",
        "stringr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  # TODO: generalize create_deg_tables?
  targets::tar_target(
    name = down_tables,
    command = create_deg_tables(
      deg_res           = res,
      comparison_list   = processed_data[["comparisons"]],
      grouping_variable = project_params[["comparison_grouping_variable"]],
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

  # This is seemly stupid as all hell, but
  # every time I attempt to get targets to iterate over
  # res itself (with `map(res)`), it converts it into a
  # list, which of course fails.  So, instead, we will
  # iterate over a vector of indices
  targets::tar_target(
    name = res_length,
    command = seq_along(res)
  ),

  targets::tar_target(
    name = unnamed_up_enrichment,
    command = deg_pathway_enrichment(
      results     = res[[res_length]],
      fcThreshold = project_params[["deg_substantial_threshold"]],
      fcPvalue    = project_params[["deg_pval_threshold"]],
      species     = project_params[["sample_species"]],
      direction   = "up",
      ontology    = "ALL"
    ),
    pattern   = map(res_length),
    iteration = "list",
    packages =
      c(
        "clusterProfiler",
        "dplyr",
        "ReactomePA"
      )
  ),

  targets::tar_target(
    name = up_enrichment,
    command = rlang::set_names(x = unnamed_up_enrichment, nm = names(res))
  ),

  targets::tar_target(
    name = unnamed_up_enrichment_degs,
    command =
      get_enrichment_fcs(
        enrichResult = up_enrichment[[res_length]],
        degResult = res[[res_length]]
          ),
    pattern = map(res_length),
    iteration = "list",
    packages = "dplyr"
  ),

  targets::tar_target(
    name = up_enrichment_degs,
    command = rlang::set_names(x = unnamed_up_enrichment_degs, nm = names(res))
  ),

  targets::tar_target(
    name = unnamed_down_enrichment,
    command = deg_pathway_enrichment(
      results     = res[[res_length]],
      fcThreshold = project_params[["deg_substantial_threshold"]],
      fcPvalue    = project_params[["deg_pval_threshold"]],
      species     = project_params[["sample_species"]],
      direction   = "down",
      ontology    = "ALL"
    ),
    pattern   = map(res_length),
    iteration = "list",
    packages =
      c(
        "clusterProfiler",
        "dplyr",
        "ReactomePA"
      )
  ),

  targets::tar_target(
    name = down_enrichment,
    command = rlang::set_names(x = unnamed_down_enrichment, nm = names(res))
  ),

  targets::tar_target(
    name = unnamed_down_enrichment_degs,
    command =
      get_enrichment_fcs(
        enrichResult = down_enrichment[[res_length]],
        degResult = res[[res_length]]
      ),
    pattern = map(res_length),
    iteration = "list",
    packages = "dplyr"
  ),

  targets::tar_target(
    name = down_enrichment_degs,
    command = rlang::set_names(x = unnamed_down_enrichment_degs, nm = names(res))
  ),

  targets::tar_target(
    name = up_tables,
    command  = create_deg_tables(
      deg_res           = res,
      comparison_list   = processed_data[["comparisons"]],
      grouping_variable = project_params[["comparison_grouping_variable"]],
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

  targets::tar_target(
    name = degs,
    command  = extract_de_genes(
      results           = res,
      comparison_list   = processed_data[["comparisons"]],
      grouping_variable = project_params[["comparison_groups"]]
    ),
    packages =
      c(
        "purrr",
        "tibble",
        "dplyr",
        "rlang"
      )
  ),

  targets::tar_target(
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

  targets::tar_target(
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

  targets::tar_target(
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

  targets::tar_target(
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

  targets::tar_target(
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

  targets::tar_target(
    name = wgcna_modules,
    command = WGCNA::blockwiseModules(
      datExpr           = vsc_top,
      power             = find_softPower(sft),
      maxBlockSize      = 20000,
      mergeCutHeight    = 0.3,
      minModuleSize     = 20,
      pamStage          = TRUE,
      pamRespectsDendro = TRUE,
      saveTOMs          = FALSE,
      verbose           = 3,
      detectCutHeight   = 0.99,
      TOMDenom          = "min",
      networkType       = "signed hybrid",
      reassignThreshold = 1e-6,
      corType = "bicor"
    ),
    packages = "WGCNA",
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = wgcna_module_genes,
    command  =
      tibble::enframe(
        x     = wgcna_modules[["colors"]],
        name = "gene",
        value = "module"
      ),
    packages = "tibble"
  ),

  targets::tar_target(
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

  targets::tar_target(
    name = wgcna_hub_genes,
    command =
      WGCNA::chooseTopHubInEachModule(
        datExpr = vsc_top,
        colorh  = wgcna_modules[["colors"]],
        power   = 4,
        type    = "signed hybrid"
      ),
    packages = "WGCNA"
  ),

  targets::tar_target(
    name = wgcna_scores,
    command  =
      dplyr::select(
        .data =
          tibble::as_tibble(
            x        = wgcna_modules[["MEs"]],
            rownames = "sample_name"
          ),
        -MEgrey
      ),
    packages =
      c(
        "dplyr",
        "tibble"
      )
  ),

  targets::tar_target(
    name = wgcna_scores_with_md,
    command =
      dplyr::left_join(
        x = wgcna_scores,
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

  targets::tar_target(
    name = unnamed_wgcna_rf_models,
    command   =
      rf_classifier(
        dataset            = wgcna_scores_with_md,
        classification_var = comparison_groups,
        starts_with("ME"),
        train_proportion   = 0.75,
        print              = FALSE
      ),
    pattern   = map(comparison_groups),
    iteration = "list",
    packages  =
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
      )
  ),

  targets::tar_target(
    name = wgcna_rf_models,
    command = rlang::set_names(x = unnamed_wgcna_rf_models, nm = comparison_groups),
    packages = "rlang"
  ),

  targets::tar_target(
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

  targets::tar_target(
    name = MEenriched_list,
    command = module_gsea(
      module_genes       = filtered_wgcna_module_genes,
      module_of_interest = wgcna_module_colors,
      target_species     = sample_species_org
    ),
    pattern   = map(wgcna_module_colors),
    iteration = "list",
    packages  =
      c(
        "dplyr",
        "clusterProfiler"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = MEenriched,
    command  = dplyr::bind_rows(MEenriched_list),
    packages = "dplyr"
  ),

  targets::tar_target(
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

  targets::tar_target(
    name = md_with_module_scores,
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

  targets::tar_target(
    name = unnamed_module_rf_models,
    command   =
      {
        message(comparison_groups)
        rf_classifier(
          dataset            = md_with_module_scores,
          classification_var = comparison_groups,
          tidyselect::matches("^M[[:digit:]]+"),
          train_proportion   = 0.75,
          print              = FALSE
        )
        },
    pattern   = map(comparison_groups),
    iteration = "list",
    packages  =
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
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = module_rf_models,
    command = rlang::set_names(x = unnamed_module_rf_models, nm = comparison_groups)
  ),

  targets::tar_target(
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

  targets::tar_target(
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

  targets::tar_target(
    name = all_module_scores,
    command =
      purrr::reduce(
        .x =
          list(
            module_scores,
            wgcna_scores
          ),
        .f = dplyr::left_join
      ),
    packages =
      c(
        "dplyr",
        "purrr"
      )
  ),

  targets::tar_target(
    name = all_module_scores_with_md,
    command =
      dplyr::left_join(
        study_md,
        all_module_scores
      ),
    packages = "dplyr"
  ),

  targets::tar_target(
    name = renamed_annotated_module_scores,
    command =
      rlang::set_names(
        nm =
          targets_recode(
            target_list = names(all_module_scores_with_md),
            thing_to_unquote_splice = annotated_mod_list
            ),
        x  = all_module_scores_with_md
      ),
    packages =
      c(
        "rlang",
        "dplyr"
      )
  ),

  targets::tar_target(
    name = module_scores_mat,
    command =
      module_scores_mat <-
      all_module_scores |>
      dplyr::select(
        sample_name,
        tidyselect::starts_with("ME"),
        tidyselect::one_of(
          c(
            names(ldg_modules),
            names(banchereau_modules),
            names(metasignature_module)
          )
        )
      ) |>
      tibble::column_to_rownames("sample_name") |>
      as.matrix() |>
      t(),
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tibble"
      )
  ),

  targets::tar_target(
    name = module_scores_pivot,
    command =
      dplyr::select(
        .data = all_module_scores_with_md,
        sample_name,
        tidyselect::matches("^M[[:digit:]]+"),
        tidyselect::matches("mg"),
        tidyselect::starts_with("ldg"),
        tidyselect::one_of(comparison_groups)
      ) %>%
      tidyr::pivot_longer(
        cols          = c(
          tidyselect::matches("^M[[:digit:]]+"),
          tidyselect::matches("mg"),
          tidyselect::starts_with("ldg")
        ),
        names_to      = "module",
        values_to     = "score"
      ),
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tidyr"
      )
  ),

  targets::tar_target(
    name = annotated_module_scores_pivot,
    command =
      dplyr::filter(
        .data = module_scores_pivot,
        module %in% annotated_modules[["module"]]
        ),
    packages = "dplyr"
  ),

  targets::tar_target(
    name = wgcna_scores_pivot,
    command =
      dplyr::select(
        .data = all_module_scores_with_md,
        sample_name,
        tidyselect::starts_with("ME"),
        tidyselect::one_of(comparison_groups)
      ) %>%
      tidyr::pivot_longer(
        cols          = c(
          tidyselect::starts_with("ME")
        ),
        names_to      = "module",
        values_to     = "score"
      ),
    packages =
      c(
        "dplyr",
        "tidyselect",
        "tidyr"
      )
  ),

  targets::tar_target(
    name = module_comparisons_stats_unnamed,
    command =
      modules_compare_with_stats(
        module_score_table = module_scores_pivot,
        compare_by = comparison_groups
      ),
    pattern = map(comparison_groups),
    iteration = "list",
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = module_comparisons_stats,
    command =
      rlang::set_names(
        x = module_comparisons_stats_unnamed,
        nm = comparison_groups
      ),
    packages = "rlang"
  ),

  targets::tar_target(
    name = wgcna_comparisons_stats_unnamed,
    command =
      modules_compare_with_stats(
        module_score_table = wgcna_scores_pivot,
        compare_by = comparison_groups
      ),
    pattern = map(comparison_groups),
    iteration = "list",
    packages =
      c(
        "dplyr",
        "stringr",
        "rstatix",
        "purrr"
      ),
    cue = targets::tar_cue(mode = "never")
  ),

  targets::tar_target(
    name = wgcna_comparisons_stats,
    command =
      rlang::set_names(
        x = wgcna_comparisons_stats_unnamed,
        nm = comparison_groups
      ),
    packages = "rlang"
  ),

  targets::tar_target(
    name = annotated_module_comparisons_stats,
    command =
      purrr::map(
        .x = module_comparisons_stats,
        .f = dplyr::filter,
        module %in% annotated_modules$module
        ),
    packages =
      c(
        "dplyr",
        "purrr"
      )
  ),

  targets::tar_target(
    name = group_pal,
    command =
      c(
        generatePalettes(
          .data = study_md,
          c(
            project_params[["comparison_grouping_variable"]],
            one_of(project_params[["heatmap_row_annotations"]])
          )
        ),
        generatePalettes(
          .data = module_annotation,
          .cols = "type",
          use_palettes = "ggsci::default_ucscgb"
          ),
        list(
          "chr" = paletteer::paletteer_d(
            palette = getRandomPalette(2),
            n = 2
            ) |>
          as.character() |>
          rlang::set_names(nm = c("X", "Y"))
          )
        ),
    packages =
      c(
        "tidyselect",
        "rlang",
        "dplyr",
        "tidyr",
        "purrr",
        "paletteer",
        "tibble"
      )
  ),

  targets::tar_target(
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

  targets::tar_target(
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

  targets::tar_target(
    name = output_module_scores,
    command  =
      save_table_to_disk(
        file_to_output = module_scores,
        output_name = "processed_data/module_scores.csv.gz"
      ),
    packages = "data.table",
    format   = "file"
  ),

  targets::tar_target(
    name = output_wgcna_scores,
    command  =
      writeData(
        object = wgcna_scores,
        output_name = "processed_data/wgcna_scores.csv.gz"
      ),
    packages = "data.table",
    format   = "file"
  ),

  # targets::tar_target(
  #   name = output_pathway_eigenvalues,
  #   command  =
  #     save_table_to_disk(
  #       file_to_output = pathway_eigenvalues,
  #       output_name = "processed_data/plafker_pathway_scores.csv.gz"
  #     ),
  #   format   = "file"
  # ),

  targets::tar_target(
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
      set_author = "Miles Smith",
      set_title  = "Shakedown run"
    ),
    output_dir    = "reports/",
    packages      =
      c(
        "rmarkdown",
        "knitr",
        "targets",
        "here",
        "ggplot2",
        "ggpubr",
        "rlang",
        "magrittr",
        "janitor",
        "kableExtra",
        "flextable",
        "genefilter",
        "tidyselect",
        "purrr",
        "formattable",
        "grid"
    ),
    quiet = FALSE
  )

  # tar_render(
  #   name = test_report,
  #   path          = "analysis/testing.rmd",
  #   params        = list(
  #     set_author = "Miles Smith",
  #     set_title  = "Shakedown run"
  #   ),
  #   output_dir    = "reports/",
  #   packages      =
  #     c(
  #       "rmarkdown",
  #       "knitr",
  #       "targets",
  #       "here",
  #       "ggplot2",
  #       "ggpubr",
  #       "rlang",
  #       "magrittr",
  #       "janitor",
  #       "kableExtra",
  #       "flextable",
  #       "genefilter",
  #       "tidyselect",
  #       "purrr",
  #       "formattable"
  #     ),
  #   quiet = FALSE
  # )
#
#   tar_render(
#     name = qc_report,
#     path          =  "analysis/qc_report.rmd",
#     params        = list(
#       set_title   = "BLAST Optimal Responder-vs-Non-responder RNAseq QC",
#       set_author  = "Miles Smith"
#     ),
#     output_dir    = "reports/"
#   )

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
