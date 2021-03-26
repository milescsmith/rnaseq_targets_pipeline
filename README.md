# Description

An R/Drake pipeline focused on analysis of RNA-sequencing.  This pipeline expects the input data to be pseudocounts from Salmon.  For a pipeline capable of processing of raw data from bcl up through alignment and quantification, see either (milescsmith/rnaseq_pipeline)[https://github.com/milescsmith/rnaseq_pipeline] for a Snakemake version or (milescsmith/nf-rnaseq)[https://github.com/milescsmith/nf-rnaseq] for a Nextflow version.

## Layout

Files in the project directory are sorted into the following subdirectories:
- `data`: Normalized data, module scores, DEG results produced during the analysis
- `external data`: Gene ontology data, transcript reference, and gene module descriptions from outside sources
- `markdown`: Rmarkdown file used to generate final report
- `metadata`: Descriptions of the samples and clinical data
- `R`: Code used for the analysis
-- `functions.R`: Custom functions used during the analysis
-- `packages.R`: All of the external libraries that need to be loaded
-- `plan.R`: the [Drake](https://github.com/ropensci/drake) analysis plan.  This describes a directed acyclic graph of *targets* that need to be produced from the input data, with the order they are produced being inferred from their interdependence (i.e. if target C requires target B to be produced, the code to produce target B will be ran first).  These targets are placed in a storr directory (the hidden `.drake` directory) as blobs that are loaded by the Rmarkdown analysis file.  Subsequent runs of the pipeline will examine the DAG and determine which portions are out of date and run only the code to produce those targets.
- `.`: contains the `_drake.R` file, which is executed by `drake::r_make()`.  Is used to source the `functions.R` and `packages.R` files, setup any global variables, and run the plan.


## Goals

The pipeline covers:
- Differential gene analysis of ? vs Control subjects
- Score previously described gene modules (Banchereau, et al. 2016 10.1016/j.cell.2016.03.008; Kegerreis, et al. 2019 10.4049/jimmunol.1801512; Haynes, et al. 2019  10.1101/834093)
- Identify patient clusters
- Identify gene expression modules associated with particular clusters, disease-classifications


## Usage

Modify the `_drake.r` file so that variables in the `#### Setup file locations ####` point to the correct locations.
In the directory containing this project, run:
  ```
docker run -it --rm --mount type=bind,source=<location_of_files>,target=<location_entered_in_drake.r> us-central1-docker.pkg.dev/scrna-196615/utopia-planitia/workerbee-rnaseq:4.0.4 /bin/Rscript "drake::r_make()"
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
