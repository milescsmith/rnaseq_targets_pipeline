sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq"
metadata_file                   = "metadata/COVID (PCV, OSCTR) analysis dataset.xlsx"
metadata_sheet                  = "main"
main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx"
main_sample_sheet               = "main"

#### Import code and libraries ####

# source("code/packages.R")        # Load your packages, e.g. library(drake).
source("code/functions.R")       # Define your custom code as a bunch of functions.
# source("code/extant_modules.R")  # Banchereau, Kegerreis, and Metagene modules


c5 <- clusterProfiler::read.gmt("references/c5.all.v6.2.symbols.gmt")

#### Set options ####
options(future.globals.maxSize = +Inf)
Sys.setenv('RSTUDIO_PANDOC'     = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()

BPPARAM                         = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)

#### Setup project variables ####
project_groups_to_include       = c("control", "OSCTR Case")
project_groups_to_exclude       = "PCV Case"
disease_classes_to_include      = NULL
disease_classes_to_exclude      = NULL
study_design                    = ~ disease_class
comparison_grouping_variable    = "disease_class"
batch_variable                  = "run_id"
control_group                   = "control"
experimental_group              = "infected"
num_sva                         = 3
comparison_groups               = c("infection_status", "hospitalization_status")

initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2
