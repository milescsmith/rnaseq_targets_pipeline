options(future.globals.maxSize = +Inf)
reticulate::use_condaenv('reticulate', required = TRUE, conda = "~/conda/bin/conda")
source("R/packages.R")  # Load your packages, e.g. library(drake).
source("R/functions.R") # Define your custom code as a bunch of functions.
c5 <- read.gmt("external_data/c5.all.v6.2.symbols.gmt")
Sys.setenv('RSTUDIO_PANDOC' = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()
### Setup bucket access
flyio_set_datasource("gcs")
flyio_auth("/home/rstudio/google_storage_access_key_scrna-196615.json")
flyio_set_bucket("memory-beta", data_source="gcs")

import_rda(file="references/gencode_v32_virus_tx2gene_v1.2.RData",
           bucket = "memory-beta") #generated from rtracklayer::readGFF()

import_rda(file="references/banchereau-modules.RData",
           bucket = "memory-beta",
           data_source = "gcs")

import_rda(file="references/kegerreis-ldg-modules.RData",
           bucket = "memory-beta",
           data_source = "gcs")

### Setup project variables
projects_to_include = c("General","MS")
projects_to_exclude = c("ALE06", "Xencor", "BChong2019.1")
disease_classes_to_include = c("Control", "MS", "SPMS", "RRMS", "NMO", "PPMS", "MS-like")
disease_classes_to_exclude = NULL
study_design = ~ sex + disease_grouping
comparison_grouping_variable = "disease_grouping"
control_group = "Control"
experimental_group = "MS_group"

initial_concentration_threshold = 1.5
pc1_zscore_threshold = 2
pc2_zscore_threshold = 2.5

### Setup file locations
seq_file_directory = "/home/rstudio/workspace/datasets/rnaseq/novaseq"
metadata_file = "metadata/NovaSeq_Sample_List.xlsx"

BPPARAM = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(future::multisession)

source("R/plan.R")
# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(plan = analysis_plan,
             verbose = 2,
             parallelism = "future",
             jobs = parallel::detectCores(),
             lock_envir = FALSE)
