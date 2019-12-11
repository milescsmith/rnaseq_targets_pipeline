options(future.globals.maxSize = +Inf)
source("R/packages.R")  # Load your packages, e.g. library(drake).
source("R/functions.R") # Define your custom code as a bunch of functions.
Sys.setenv('RSTUDIO_PANDOC' = '/usr/lib/rstudio-server/bin/pandoc')
### Setup bucket access
flyio_set_datasource("gcs")
flyio_auth("/opt/google_project_scrna_196615_storage_key.json")
flyio_set_bucket("memory-beta", data_source="gcs")

import_rda(file="references/gencode_v32_virus_tx2gene.RData",
           bucket = "memory-beta") #generated from rtracklayer::readGFF()

import_rda(file="references/banchereau-modules.RData",
           bucket = "memory-beta",
           data_source = "gcs")

import_rda(file="references/kegerreis-ldg-modules.RData",
           bucket = "memory-beta",
           data_source = "gcs")

### Setup project variables
projects_to_include = c("General","DxTerity","ABC","MS","OCRD","MDFL","ORDRCC.13-19","ORDRCC", "BLAST", "OMRF.13-35")
projects_to_exclude = c("ALE06", "Xencor", "BChong2019.1")
disease_classes_to_include = c("Control", "SLE")
disease_classes_to_exclude = NULL
study_design = ~ initial_concentration_ng_ul + run_id + disease_class
comparison_grouping_variable = "disease_class"
control_group = "Control"
experimental_group = "SLE"

initial_concentration_threshold = 1.5
pc1_zscore_threshold = 2

### Setup file locations
seq_file_directory = "/media/charon/datasets/viral_transcripts"
metadata_file = "datasets/rnaseq/novaseq/NovaSeq_Sample_List.xlsx"

BPPARAM = SnowParam(workers=parallel::detectCores()-4, type = "SOCK")
register(BPPARAM)

source("R/plan.R")
# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(plan,
             verbose = 2,
             parallelism = "future",
             jobs = parallel::detectCores()-4,
             lock_envir = FALSE)