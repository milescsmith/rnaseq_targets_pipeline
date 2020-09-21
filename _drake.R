options(future.globals.maxSize = +Inf)
source("R/packages.R")  # Load your packages, e.g. library(drake).
source("R/functions.R") # Define your custom code as a bunch of functions.
c5 <- read.gmt("external_data/c5.all.v6.2.symbols.gmt")
Sys.setenv('RSTUDIO_PANDOC' = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()

load(file="external_data/gencode_v32_virus_tx2gene_v1.2.RData") #generated from rtracklayer::readGFF()

load(file="external_data/banchereau-modules.RData")

load(file="external_data/kegerreis-ldg-modules.RData")

metasignature_module <- list(mg = c("IFI44L", "EPSTI1", "HERC5",
                        "IFI44", "IFI27", "RSAD2",
                        "ISG15", "OASL", "IFIT3",
                        "CMPK2", "IFIT1", "USP18",
                        "SIGLEC1", "MX2", "MX1",
                        "LY6E", "PLSCR1", "OAS2",
                        "SPATS2L", "OAS1", "OAS3",
                        "PARP12", "SCO2", "IFIH1",
                        "IFITM3", "IFITM1", "DDX60",
                        "TRIM22", "RTP4", "SAMD9L",
                        "XAF1", "IFI35", "TNFAIP6",
                        "MT2A", "LAP3", "HERC6",
                        "FBXO6", "TDRD7", "IFIT5",
                        "PHF11", "IFIT2", "TAP1",
                        "TOR1B", "TNFSF13B", "ELANE",
                        "IRF7", "DHX58", "ZBP1",
                        "TYMP", "GRN", "LAMP3",
                        "IFI6", "NTNG2", "SAT1",
                        "SERPING1", "DEFA4", "DDX58",
                        "SAMD9", "PARP9", "MT1E",
                        "IRF9", "TCN2", "MT1HL1",
                        "HSH2D", "ISG20", "ZCCHC2",
                        "SP100", "LMO2", "GBP1",
                        "IFITM2", "LGALS3BP", "REM2",
                        "TNFSF10", "HESX1", "MT1A",
                        "CCR1", "CEACAM1", "MYD88",
                        "BST2", "LHFPL2", "FFAR2",
                        "MT1F", "RPL10A", "CD1C",
                        "VSIG1", "GPR183", "DSC1",
                        "KLRB1", "FBL", "EIF3E",
                        "ABCB1", "NAP1L3", "EIF3L"))

cytokine_names <- c("baff" = "BAFF",
                    "blc_cxcl13" = "BCL-CXCL13",
                    "gm_csf" = "GM-CSF",
                    "gro_alpha" = "GROalpha",
                    "ifn_alpha" = "IFNalpha",
                    "ifn_gamma" = "IFNgamma",
                    "il1alpha" = "IL1alpha",
                    "il1beta" = "IL1beta",
                    "il_1ra" = "IL1RA",
                    "il_10" = "IL10",
                    "il_12p70" = "IL12p70",
                    "il_13" = "IL13",
                    "il_15" = "IL15",
                    "il_17a" = "IL17a",
                    "il_18" = "IL18",
                    "il_2" = "IL2",
                    "il_21" = "IL21",
                    "il22" = "IL22",
                    "il_23" = "IL23",
                    "il_27" = "IL27",
                    "il_31" = "IL31",
                    "il_4" = "IL4",
                    "il5" = "IL5",
                    "il_7" = "IL7",
                    "il_8" = "IL8",
                    "il_9" = "IL9",
                    "ip_10" = "IP10",
                    "mcp1" = "MCP1",
                    "mip1alpha" = "MIP1alpha",
                    "mip_1beta" = "MIP1beta",
                    "rantes" = "RANTES",
                    "s_cd40l" = "sCD40L",
                    "scf" = "SCF",
                    "sdf_1alpha" = "SDF1alpha",
                    "s_icam_1" = "sICAM1",
                    "tnf_alpha" = "TNFalpha",
                    "tnf_beta" = "TNFbeta")

### Setup project variables
projects_to_include = "BChong2019.1"
projects_to_exclude = c("ALE06", "Xencor")
disease_classes_to_include = NULL
disease_classes_to_exclude = NULL
samples_to_remove = c("UTSW-BC02","UTSW-BC29","UTSW-BC48", "UTSW_BC38")
study_design = ~ disease_class
comparison_grouping_variable = "disease_class"
control_group = "Control"
experimental_group = "DLE"
batch_variable = "run_id"
num_sva = 3

initial_concentration_threshold = 1.5
pc1_zscore_threshold = 2
pc2_zscore_threshold = 2.5

### Setup file locations
seq_file_directory = "/home/rstudio/workspace/datasets/rnaseq/novaseq"
metadata_file = "metadata/NovaSeq_Sample_List.xlsx"
clinical_file = "metadata/study_data.csv"
sledai_bpx_file = "metadata/SLEDAI_BPX.xlsx"

BPPARAM = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(future::multisession)

source("R/plan.R")
# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(plan = analysis_plan,
             verbose = 2,
             parallelism = "future",
             jobs = parallel::detectCores()-2,
             lock_envir = FALSE)
