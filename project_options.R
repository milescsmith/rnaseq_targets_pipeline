seq_file_directory              = "/home/rstudio/oldworkspace/datasets/rnaseq/novaseq"
metadata_file                   = "metadata/COVID (PCV, OSCTR) analysis dataset.xlsx"
metadata_sheet                  = "main"
main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx"
main_sample_sheet               = "main"

#### Import code and libraries ####

source("code/packages.R")        # Load your packages, e.g. library(drake).
source("code/functions.R")       # Define your custom code as a bunch of functions.
source("code/extant_modules.R")  # Banchereau, Kegerreis, and Metagene modules


c5 <- read.gmt("references/c5.all.v6.2.symbols.gmt")

#### Set options ####
options(future.globals.maxSize = +Inf)
Sys.setenv('RSTUDIO_PANDOC'     = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()


#### Setup project variables ####
project_groups_to_include       = c("PCV Case", "control")
project_groups_to_exclude       = c("OSCTR Case")
disease_classes_to_include      = NULL
disease_classes_to_exclude      = NULL
study_design                    = ~ run_id + disease_class
comparison_grouping_variable    = "disease_class"
batch_variable                  = "run_id"
control_group                   = "control"
experimental_group              = "infected"
num_sva                         = 3

initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2
