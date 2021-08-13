sequencing_file_directory       = "/home/rstudio/workspace/datasets/rnaseq/novaseq"
metadata_file                   = "metadata/old/blast_metadata.csv"
main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx"
main_sample_sheet               = "main"

#### Import code and libraries ####
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

study_design                    = ~ responder
comparison_grouping_variable    = "responder"
batch_variable                  = "run_id"
control_group                   = "non_responder"
experimental_group              = "responder"
num_sva                         = 3
timepoint                       = "BL"
manual_sample_removal           = NULL #c("AA04629", "AA04874")
responders                      = c(1, 2, 5,13,17,34)
non_responders                  = c(7,10,23,24,26,31)
groups_to_include               = c("Unaffected Control", "Vaccinated", "Infected + Vaccine", "Infected")
groups_to_exclude               = c("OSCTR Case")
initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2.5

# number of top variable genes to use for WGCNA
n_var_genes                     = 20000
