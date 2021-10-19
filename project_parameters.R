#### Set options ####
options(future.globals.maxSize    = +Inf)
Sys.setenv('RSTUDIO_PANDOC'       = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()
BPPARAM =
  BiocParallel::SnowParam(
    workers       = parallel::detectCores()-2,
    exportglobals = FALSE,
    progressbar   = TRUE,
    type = "SOCK"
    )

BiocParallel::register(BPPARAM)


#### Setup project variables ####
project_params = list(

  sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq",
  metadata_file                   = "metadata/blast_baseline.csv",
  metadata_sheet                  = NULL,
  skip_lines                      = 0,
  main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx",
  main_sample_sheet               = NULL,
  main_sample_sheet_skip          = 1,
  annotation_file                 = "references/gencode_v32_virus_tx2gene_v1.2.csv",

  #### Metadata columns ####
  sample_species                  = "Homo sapiens",
  sample_name_column              = "NovaSeq_Sample_ID",
  grouping_column                 = "responder",
  project_column                  = "project",
  regression_columns              = c("age", "sex", "RaceCode", "final_concentration_ng_ul"),
  filter_column                   = "initial_concentration_ng_ul",
  filter_value                    = 1.5,
  extra_columns                   = NULL,
  heatmap_row_annotations         = c("responder", "sex",  "RaceCode", "cluster"),

  #### Setup project variables ####
  projects_to_include             = "BLAST",
  projects_to_exclude             = c("ALE06", "Xencor", "BChong2019.1"),

  groups_to_include               = c("responder", "non_responder"),
  groups_to_exclude               = NULL,
  study_design                    = ~ responder,

  comparison_grouping_variable    = "responder",
  comparison_groups               = "responder",
  batch_variable                  = "run_id",
  control_group                   = "non_responder",
  experimental_group              = "responder",
  manual_sample_removal           = c("AA04629", "AA04874"),

  aligner                         = "salmon",
  only_hugo_named_genes           = FALSE,
  minimum_gene_count              = 1,
  initial_concentration_threshold = 1.5,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2.5,
  sva_num                         = 3,
  use_combat                      = FALSE,
  process_method                  = "DESeq2",
  lfcThreshold                    = 0.01,
  deg_substantial_threshold       = 0.01, # which is about a 2.5-fold change
  deg_pval_threshold              = 0.05,

  BPPARAM                         = BPPARAM,

  # number of top variable genes to use for WGCNA
  n_var_genes                     = 20000,

  banchereau_modules              = "references/banchereau_modules.csv",
  banchereau_module_annotations   = "references/banchereau_module_annotations.csv",
  ldg_modules                     = "references/ldg_modules.csv",
  metasignature_modules           = "references/metasignature_module.csv"
  )
