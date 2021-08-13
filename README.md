{[targets](https://docs.ropensci.org/targets/)}-based RNAseq analysis pipeline

# Warning

Some aspects of this still have hardcoded variable names that are a hold over from my previous {[drake](https://github.com/ropensci/drake)}-based [rnaseq_drake](https://github.com/milescsmith/rnaseq_drake) pipeline.  Consider this still something of a work-in-progress.

# Metadata

You will need to supply a metadata file describing the samples for processing.
At a minimum, you will need the following columns:
* sample_name - this should match (at least part) of the quant.sf.gz filename or path to it.
* sex
* ethnicity
* age
* project
* study_group
* rin
* concentration


## Goals

The pipeline covers:
- Differential gene analysis of a defined experimental group vs Control subjects
- Score previously described gene modules (Banchereau, et al. 2016 10.1016/j.cell.2016.03.008; Kegerreis, et al. 2019 10.4049/jimmunol.1801512; Haynes, et al. 2019  10.1101/834093)
- Identify sample clusters using random forest proximities and k-means clustering
- Identify gene expression modules associated with particular clusters, disease-classifications

## Structure

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

## Usage

Modify the `project_options.r` file so that variables in the `### Setup file locations` point to the correct locations.
In the directory containing this project, run:
  ```
docker run -it --rm --mount type=bind,source=<location_of_files>,target=<location_entered_in_drake.r> us-central1-docker.pkg.dev/guthridge-nih-strides-projects/utopia-planitia/workerbee-rnaseq:4.1.0b /bin/Rscript "targets::tar_make()"
```

(Note: this analysis may require a lot of memory.  It has been sucessfully run with 64G, but I am unclear how much is actually required.)

### Analysis on the OMRF HPC

Since Docker is incompatible with use on a shared compute cluster (Docker requires elevated permissions that no sane admin would grant), it is necessary to use [Singularity](https://sylabs.io/guides/3.5/user-guide/).  To run:
  ```
module load singularity
singularity exec --bind=<location_entered_in_drake.r> /Volumes/guth_aci_informatics/software/workerbee-rnaseq-drake-1.6.sif Rscript -e "drake::r_make('<path to _drake.r>')"
```

## Last Updated
See the [CHANGELOG](CHANGELOG.md) for the latest updates
