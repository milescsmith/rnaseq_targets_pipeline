sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/narch_advanta/data/results/aligned/"
metadata_file                   = "metadata/narch_rnaseq_targets_metadata.csv"

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
projects_to_include             = NULL
projects_to_exclude             = NULL
study_groups_to_include         = c("control", "sle", "ra")
study_groups_to_exclude         = NULL
study_design                    = ~ plate_id + study_group
comparison_grouping_variable    = "study_group"
batch_variable                  = "plate_id"
control_group                   = "control"
experimental_group              = c("sle", "ra")
num_sva                         = 3
comparison_groups               = c("study_group", "sample_type")

initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2
