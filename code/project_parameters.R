#### Set options ####
options(future.globals.maxSize    = +Inf)
Sys.setenv('RSTUDIO_PANDOC'       = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()
BPPARAM =
  BiocParallel::MulticoreParam(
    workers       = parallel::detectCores()-1,
    exportglobals = FALSE,
    progressbar   = TRUE
    )

BiocParallel::register(BPPARAM)


#### Setup project variables ####
project_params = list(

  sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq",
  metadata_file                   = "metadata/blast_metadata.csv",
  main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx",
  main_sample_sheet               = "main",
  annotation_file                 = "references/gencode_v32_virus_tx2gene_v1.2.csv",

  #### Metadata columns ####
  sample_species                  = "Homo sapiens",
  sample_name_column              = "NovaSeq_Sample_ID",
  grouping_column                 = "Disease_Class",
  project_column                  = "project",
  regression_columns              = c("age", "sex", "run_id", "RaceCode"),
  filter_column                   = "initial_concentration_ng_ul",
  filter_value                    = 1.5,
  heatmap_row_annotations         = c("Disease_Class", "sex", "run_id", "RaceCode", "cluster"),

  #### Setup project variables ####
  projects_to_include             = "BLAST",
  projects_to_exclude             = c("Xencor", "DxTerity"),

  groups_to_include               = c("Control", "SLE"),
  groups_to_exclude               = NULL,
  study_design                    = NULL, # ~ Disease_Class,

  comparison_grouping_variable    = "Disease_Class",
  comparison_groups               = c("Disease_Class", "cluster"),
  batch_variable                  = "run_id",
  control_group                   = "control",
  experimental_group              = "sle",
  num_sva                         = 3,
  manual_sample_removal           = NULL,

  aligner                         = "salmon",
  only_hugo_named_genes           = TRUE,
  minimum_gene_count              = 1,
  initial_concentration_threshold = 1.5,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2.5,
  sva_num                         = 2,
  use_combat                      = FALSE,
  process_method                  = "limma",
  lfcThreshold                    = 0.25,

  BPPARAM                         = BPPARAM,

  # number of top variable genes to use for WGCNA
  n_var_genes                     = 20000,

  banchereau_modules              = "references/banchereau_modules.csv",
  banchereau_module_annotations   = "references/banchereau_module_annotations.csv",
  ldg_modules                     = "references/ldg_modules.csv",
  metasignature_modules           = "references/metasignature_module.csv",

  # Report parameters
  row_annotations                 = c("Disease_Class", "project", "sex", "run_id", "RaceCode", "cluster")
  )
