seq_file_directory              = "/home/rstudio/oldworkspace/datasets/rnaseq/novaseq"
metadata_file                   = "metadata/COVID (PCV, OSCTR) analysis dataset.xlsx"
metadata_sheet                  = "main"
main_sample_list                = "metadata/NovaSeq_Sample_List.xlsx"
main_sample_sheet               = "main"

#### Import code and libraries ####

source("code/packages.R")        # Load your packages, e.g. library(drake).
source("code/functions.R")       # Define your custom code as a bunch of functions.
source("code/extant_modules.R")  # Banchereau, Kegerreis, and Metagene modules
load(file="references/gencode_v32_virus_tx2gene_v1.2.RData") #generated from rtracklayer::readGFF()

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
study_design                    = ~ disease_class
comparison_grouping_variable    = "disease_class"
control_group                   = "control"
experimental_group              = "infected"
batch_variable                  = "run_id"
num_sva                         = 9

samples_to_manually_remove      = c("AA05181", "AA05141", "AA04447")

initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2

#### Parallelism ####
BPPARAM                         = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(future::multisession)


#### Import and assemble analysis plan parts ####
import_plan                     = code_to_plan("code/plan/01_data_import.R")
filtering_plan                  = code_to_plan("code/plan/02_filtering_and_batch_correction.R")
processing_plan                 = code_to_plan("code/plan/03_processing.R")
module_scoring_plan             = code_to_plan("code/plan/04_module_scoring.R")
clustering_plan                 = code_to_plan("code/plan/05_clustering.R")
dimensional_reduction_plan      = code_to_plan("code/plan/06_dimensional_reduction.R")
de_plan                         = code_to_plan("code/plan/07_differential_expression.R")
WGCNA_plan                      = code_to_plan("code/plan/08_WGCNA.R")
classifier_plan                 = code_to_plan("code/plan/09_module_classifier.R")
viral_plan                      = code_to_plan("code/plan/10_viral_transcripts.R")
stats_plan                      = code_to_plan("code/plan/11_stats_testing.R")
output_plan                     = code_to_plan("code/plan/12_table_outputs.R")
palettes_plan                   = code_to_plan("code/plan/98_palettes.R")
reports_plan                    = code_to_plan("code/plan/99_reports.R")

compiled_analysis_plan <-
  bind_plans(
    import_plan,
    filtering_plan,
    processing_plan,
    module_scoring_plan,
    clustering_plan,
    dimensional_reduction_plan,
    de_plan,
    WGCNA_plan,
    classifier_plan,
    viral_plan,
    stats_plan,
    output_plan,
    palettes_plan,
    reports_plan
  )

# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
#### Run plan ####
drake_config(
  plan                    = compiled_analysis_plan,
  cache_log_file          = here("logs/cache.log"),
  log_make                = here("logs/console.log"),
  log_progress            = TRUE,
  log_build_times         = TRUE,
  session_info            = TRUE,
  verbose                 = 2,
  # caching                 = "worker",
  # memory_strategy         = "speed",
  parallelism             = "future",
  jobs                    = parallel::detectCores()-2,
  lock_envir              = FALSE,
  #targets                 = c("report", "qc_report", "report_pdf", "qc_report_pdf")
  targets                 = "report"
  )
