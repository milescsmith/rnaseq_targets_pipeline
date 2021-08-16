#### Import code and libraries ####
source("code/functions.R")       # Define your custom code as a bunch of functions.

c5 <- clusterProfiler::read.gmt("references/c5.all.v7.4.symbols.gmt")

#### Set options ####
options(future.globals.maxSize    = +Inf)
Sys.setenv('RSTUDIO_PANDOC'       = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()
BPPARAM                           = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)


#### Setup project variables ####
project_params = list(

  sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq",
  metadata_file                   = "metadata/blast_metadata.csv",
  main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx",
  main_sample_sheet               = "main",
  annotation_file                 = "references/gencode_v32_virus_tx2gene_v1.2.csv",

  #### Metadata columns ####
  sample_name_column              = "NovaSeq_Sample_ID",
  grouping_column                 = "Disease_Class",
  project_column                  = "project",
  regression_columns              = c("age", "sex", "run_id", "RaceCode"),
  filter_column                   = "initial_concentration_ng_ul",
  filter_value                    = 1.5,

  #### Setup project variables ####
  projects_to_include             = c("BLAST"),
  projects_to_exclude             = NULL,

  groups_to_include               = c("Control", "SLE"),
  groups_to_exclude               = NULL,
  study_design                    = NULL, # ~ Disease_Class,

  comparison_grouping_variable    = "Disease_Class",
  batch_variable                  = "plate_id",
  control_group                   = "control",
  experimental_group              = c("sle", "ra"),
  num_sva                         = 3,
  #comparison_groups               = c("study_group", "sample_type"),
  manual_sample_removal           = NULL,

  aligner                         = "salmon",
  minimum_gene_count              = 1,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2,
  BPPARAM                         = BPPARAM,
  sva_num                         = 2,
  use_combat                      = FALSE,
  process_method                  = "DESeq2",

  initial_concentration_threshold = 1.5,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2.5,

  # number of top variable genes to use for WGCNA
  n_var_genes                     = 20000,

  # TODO: generalize the module import functions so that things are not explicit?
  banchereau_modules              = "references/banchereau_modules.csv",
  banchereau_module_annotations   = "references/banchereau_module_annotations.csv",
  ldg_modules                     = "references/ldg_modules.csv",
  metasignature_modules           = "references/metasignature_module.csv"
)
