options(future.globals.maxSize = +Inf)
source("R/packages.R")  # Load your packages, e.g. library(drake).
source("R/functions.R") # Define your custom code as a bunch of functions.
c5 <- read.gmt("references/c5.all.v6.2.symbols.gmt")
#WGCNA::allowWGCNAThreads()

load(file="references/gencode_v32_virus_tx2gene_v1.2.RData") #generated from rtracklayer::readGFF()
load(file="references/banchereau-modules.RData")
load(file="references/kegerreis-ldg-modules.RData")

### Setup project variables
projects_to_include = NULL
projects_to_exclude = c("ALE06", "Xencor", "BChong2019.1")
disease_classes_to_include = c("Control", "SLE")
disease_classes_to_exclude = NULL
study_design = ~ disease_class
comparison_grouping_variable = "disease_class"
control_group = "Control"
experimental_group = "SLE"

initial_concentration_threshold = 1.5
pc1_zscore_threshold = 2
pc2_zscore_threshold = NULL

### Setup file locations
seq_file_directory = "/home/rstudio/workspace/datasets/rnaseq/novaseq"
metadata_file = "metadata/NovaSeq_Sample_List.xlsx"

BPPARAM = BiocParallel::SnowParam(workers=48, type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(future::multisession)

source("R/plan.R")
# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(plan = analysis_plan,
             verbose = 2,
             parallelism = "future",
             log_progress = TRUE,
             jobs = 48,
             lock_envir = FALSE)
