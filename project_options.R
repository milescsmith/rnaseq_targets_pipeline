sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq"
metadata_file                   = "metadata/COVID (PCV, OSCTR) analysis dataset.xlsx"
metadata_sheet                  = "main"
main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx"
main_sample_sheet               = "main"

#### Import code and libraries ####

# source("code/packages.R")        # Load your packages, e.g. library(drake).
source("code/functions.R")       # Define your custom code as a bunch of functions.
# source("code/extant_modules.R")  # Banchereau, Kegerreis, and Metagene modules


c5 <- clusterProfiler::read.gmt("references/c5.all.v7.4.symbols.gmt")

#### Set options ####
options(future.globals.maxSize = +Inf)
Sys.setenv('RSTUDIO_PANDOC'     = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()

BPPARAM                         = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)

#### Setup project variables ####
groups_to_include               = c("Unaffected Control", "Vaccinated", "Infected + Vaccine", "Infected")
groups_to_exclude               = c("OSCTR Case")
study_design                    = ~ study_group
comparison_grouping_variable    = "study_group"
batch_variable                  = "run_id"
control_group                   = "Unaffected Control"
# experimental_group              = "infected"
num_sva                         = 3
manual_sample_removal           = c("AA05222")

initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2
