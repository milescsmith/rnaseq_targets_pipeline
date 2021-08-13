# COVID RNAseq data - 2021/08/04

## Goals
- Examine transcriptional differences between infected, infected and vaccinated, uninfected, and uninfected and vaccinated samples
- Differential gene analysis
- Score previously described gene modules (Banchereau, et al. 2016 10.1016/j.cell.2016.03.008; Kegerreis, et al. 2019 10.4049/jimmunol.1801512; Haynes, et al. 2019  10.1101/834093)
- Identify sample clusters
- Identify gene expression modules associated with particular clusters, disease-classifications

## Analyst
- Miles Smith <miles-smith@omrf.org>

## Description
RNA-sequencing of RNA isolated from whole blood. RNA isolated by Qiagen RNeasy, libraries contructed using
Fluidigm Advanta-seq, and sequenced on a Novaseq 6000 S4.  Raw data was
converted from bcls to FASTQs using bcl2fastq 2.20, trimmed using bbmap 38.86, pseudoaligned and quantified
with Salmon 1.3.0 using a index built from GENCODE v32.  All but the bcl -> FASTQ conversion was performed
using a Nextflow script (https://gitlab.com/milothepsychic/rnaseq_pipeline/,
commit 1e2826221c220768c8cfd68fadcded43a864ea98) and docker containers.

Files in the project directory are sorted into the following subdirectories:
- `processed_data`: Normalized data, module scores, DEG results produced during the analysis
- `references`: Gene ontology data, transcript reference, and gene module descriptions from outside sources
- `analysis`: Rmarkdown file used to generate final report
- `metadata`: Descriptions of the samples and clinical data
- `code`: Code used for the analysis
-- `functions.R`: Custom functions used during the analysis
-- `packages.R`: All of the external libraries that need to be loaded
-- `plan`: Custom functions that are specific to a particular portion of the analysis pipeline.  These are broken up roughly by the stage at which they are used.
- `_targets.R`: the [targets](https://github.com/ropensci/targets) analysis plan.  This describes a directed acyclic graph of analysis *targets* that need to be produced from the input data, with the order they are produced being inferred from their interdependence (i.e. if target C requires target B to be produced, the code to produce target B will be ran first).  These targets are placed in the `_targets` directory and are loaded by the Rmarkdown analysis file.  Subsequent runs of the pipeline will examine the DAG and determine which portions are out of date and run only the code to produce those targets.

The pipeline is invoked `targets::tar_make()`.  That runs the `_targets.R` file, which in turn sources `project_options.R` and the files within the `code` directory and, by default, produces all of the non-existent or out of date targets.

## Data sources
RNA was extracted by Nicolas Dominguez and sequencing libraries
were created by Fluidigm.  Sequencing was performed in the
OMRF Clinical Genomics Core.  RNAseq data was processed by me (Miles Smith) as
detailed above.

## Usage

Modify the `project_options.r` file so that variables in the `### Setup file locations` point to the correct locations.
In the directory containing this project, run:
  ```
docker run -it --rm --mount type=bind,source=<location_of_files>,target=<location_entered_in_drake.r> us-central1-docker.pkg.dev/guthridge-nih-strides-projects/utopia-planitia/workerbee-rnaseq:4.1.0b /bin/Rscript "targets::tar_make()"
```

(Note: this analysis may require a lot of memory.  It has been successfully run with 64G, but I am unclear how much is actually required.)

### Analysis on the OMRF HPC

Since Docker is incompatible with use on a shared compute cluster (Docker requires elevated permissions that no sane admin would grant), it is necessary to use [Singularity](https://sylabs.io/guides/3.5/user-guide/).  To run:
  ```
module load singularity
singularity exec --bind=<location_entered_in_drake.r> /Volumes/guth_aci_informatics/software/workerbee-rnaseq-drake-1.6.sif Rscript -e "drake::r_make('<path to _drake.r>')"
```

## Last Updated
See the [CHANGELOG](CHANGELOG.md) for the latest updates
