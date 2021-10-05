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
  metadata_file                   = "metadata/metadata.csv",
  metadata_sheet                  = "Visits",
  skip_lines                      = 1,
  main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx",
  main_sample_sheet               = "Visits",
  main_sample_sheet_skip          = 1,
  annotation_file                 = "references/gencode_v32_virus_tx2gene_v1.2.csv",

  #### Metadata columns ####
  sample_species                  = "Homo sapiens",
  sample_name_column              = "name",
  grouping_column                 = "StudyGroup",
  project_column                  = "StudyGroup",
  regression_columns              = c("AgeAtSample", "Sex"),
  filter_column                   = NULL,
  filter_value                    = NULL,
  extra_columns                   = "RaceCode",
  heatmap_row_annotations         = c("StudyGroup", "Sex",  "RaceCode", "cluster"),

  #### Setup project variables ####
  projects_to_include             = NULL,
  projects_to_exclude             = NULL,

  groups_to_include               = c("Unaffected Control", "Infected"),
  groups_to_exclude               = NULL,
  study_design                    = ~ StudyGroup,

  comparison_grouping_variable    = "StudyGroup",
  comparison_groups               = c("StudyGroup", "cluster"),
  batch_variable                  = NULL,
  control_group                   = "control",
  experimental_group              = "sle",
  manual_sample_removal           = "#no sample",

  aligner                         = "salmon",
  only_hugo_named_genes           = TRUE,
  minimum_gene_count              = 1,
  initial_concentration_threshold = 1.5,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2.5,
  sva_num                         = 3,
  use_combat                      = FALSE,
  process_method                  = "limma",
  lfcThreshold                    = 0.25,
  deg_substantial_threshold       = 1.32, # which is about a 2.5-fold change
  deg_pval_threshold              = 0.05,

  BPPARAM                         = BPPARAM,

  # number of top variable genes to use for WGCNA
  n_var_genes                     = 20000,

  banchereau_modules              = "references/banchereau_modules.csv",
  banchereau_module_annotations   = "references/banchereau_module_annotations.csv",
  ldg_modules                     = "references/ldg_modules.csv",
  metasignature_modules           = "references/metasignature_module.csv"
  )
