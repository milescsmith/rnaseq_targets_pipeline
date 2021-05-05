#### Import code and libraries ####

source("code/functions.R")       # Define your custom code as a bunch of functions.

c5 <- clusterProfiler::read.gmt("references/c5.all.v6.2.symbols.gmt")

#### Set options ####
options(future.globals.maxSize = +Inf)
Sys.setenv('RSTUDIO_PANDOC'     = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()

BPPARAM                         = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)

#### Setup project variables ####
project_params = list(
  sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/narch_advanta/data/results/aligned/",
  metadata_file                   = "metadata/narch_rnaseq_targets_metadata.csv",
  annotation_file                 = "references/gencode_v32_virus_tx2gene_v1.2.csv",
  projects_to_include             = NULL,
  projects_to_exclude             = NULL,
  study_groups_to_include         = c("control", "sle", "ra"),
  study_groups_to_exclude         = NULL,
  study_design                    = ~ run_id + sample_type + study_group,
  comparison_grouping_variable    = "study_group",
  batch_variable                  = "run_id",
  control_group                   = "control",
  experimental_group              = c("sle", "ra"),
  num_sva                         = 3,
  comparison_groups               = c("study_group", "sample_type"), #all of these **MUST** be in the study_design

  aligner                      = "salmon",
  minimum_gene_count           = 1,
  pc1_zscore_threshold         = 2,
  pc2_zscore_threshold         = 2,
  BPPARAM                      = BPPARAM,
  sva_num                      = 2,
  use_combat                   = FALSE,
  process_method               = "DESeq2",


  initial_concentration_threshold = 1.5,
  pc1_zscore_threshold            = 2,
  pc2_zscore_threshold            = 2,

  # TODO: generalize the module import functions so that things are not explicit?
  banchereau_modules            = "references/banchereau_modules.csv",
  banchereau_module_annotations = "references/banchereau_module_annotations.csv",
  ldg_modules                   = "references/ldg_modules.csv",
  metasignature_modules         = "references/metasignature_module.csv"
)
