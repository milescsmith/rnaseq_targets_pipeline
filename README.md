{[targets](https://docs.ropensci.org/targets/)}-based RNAseq analysis pipeline

## Goals

This pipeline aims to make a generalized RNAseq analysis easier, with minimal
(or no) additional coding necessary between analyses.  Currently, it works
with sequencing files output by the pseudoaligner Salmon, but should work file
with those from Kallisto; some work will be required to enable output by STAR/
Stringtie.

The pipeline covers:
- Differential gene analysis of defined experimental group(s) vs control
subjects
-- Examination of differentially expressed genes for pathway enrichment
- Score previously described gene modules (Banchereau, et al. 2016
10.1016/j.cell.2016.03.008; Kegerreis, et al. 2019 10.4049/jimmunol.1801512; 
Haynes, et al. 2019  10.1101/834093)
- Identify sample clusters using random forest proximities and k-means
clustering
- Identify gene expression modules associated with particular clusters and/or
experimental groups using weighted gene correlation network analysis
-- Perform gene set enrichment analysis to identify pathways associated with
WGCNA modules
- Perform modeling (for now, just random forest classification) to identify
modules most important for classification

## Usage

Modify the `code/project_parameters.R` file so that variables in the 
`project_params` list
1. Point to the correct locations of the metadata and seqeucning files.
2. Identify important columns in the metadata file
3. Define experimental groups and experimental study design
4. Set the value for settings important in analysis, such as cut-off thresholds

Beyond variables present in the metadata file, probably the most important are
the pc1_zscore_threshold and pc2_zscore_threshold - initially, all desired
samples are imported and undergo a minimal amount of processing.  A PCA is then
performed on the inital normalized data, a z-score is calculated for principle
components 1 and 2, and samples with a PC1/PC2 z-score beyond the set threshold 
are removed from subsequent analysis as they are assumed to be outliers.

## Running

In the directory containing this project, run:
  ```
docker run -it --rm --mount type=bind,source=<location_of_files>,target=<values_present_in_project_parameters> us-central1-docker.pkg.dev/guthridge-nih-strides-projects/utopia-planitia/workerbee-rnaseq:4.1.1 /bin/Rscript "targets::tar_make()"
```

(Note: this analysis may require a lot of memory.  It has been successfully run 
with 64G, but I am unclear how much is actually required.)

### Metadata

The metadata file can be either a simple comma/tab-delimited file or 
and Excel spreadsheet; if it is the latter, you will need to specify the 
name of the sheet containing the metadata.  In `project_parameters.R` is a list
of the types of data required, which include:
* sample_species - either the simple name ("human") or species name
("Homo sapiens") should work.  Currently, only mice and humans are functional
* sample_name - this should match (at least part) of the sequence filename
or the path to it.
* grouping_column - should match the "comparison_grouping_variable".  Here
so that it is included in the import of metadata
* project_column - samples and sequences may be in one big pile.  This allows
for the extraction of only particular samples
* regression_columns - covariates to regress out
* filter_column - a numeric column used to remove samples of low quality (such
as a low RIN value or low RNA concentration)

### Analysis on the OMRF HPC

Since Docker is incompatible with use on a shared compute cluster (Docker requires elevated permissions that no sane admin would grant), it is necessary to use [Singularity](https://sylabs.io/guides/3.5/user-guide/).  To run:
  ```
module load singularity
singularity exec --bind=<location_entered_in_drake.r> /Volumes/guth_aci_informatics/software/workerbee-rnaseq-drake-4.1.1.sif Rscript -e "targets::tar_make()"
```

## Structure

Files in the project directory are sorted into the following subdirectories:
- `processed_data`: Normalized data, module scores, DEG results produced during the analysis
- `references`: Gene ontology data, transcript reference, and gene module descriptions from outside sources
- `analysis`: Rmarkdown file used to generate final report
-- `templates`: Child Rmarkdown files that are used in loops present in the main report.  Allow iteration over variables of unpredetermined length
- `metadata`: Descriptions of the samples and clinical data
- `code`: Code used for the analysis
-- `project_parameters.R`: a file in which parameters for the analysis should be set, eliminating the need to recode anything in `_targets.R` or in files under `plan`
-- `plan`: Custom functions that are specific to a particular portion of the analysis pipeline.  These are broken up roughly by the stage at which they are used.
- `_targets.R`: the [targets](https://github.com/ropensci/targets) analysis plan.  This describes a directed acyclic graph of analysis *targets* that need to be produced from the input data, with the order they are produced being inferred from their interdependence (i.e. if target C requires target B to be produced, the code to produce target B will be ran first).  These targets are placed in the `_targets` directory and are loaded by the Rmarkdown analysis file.  Subsequent runs of the pipeline will examine the DAG and determine which portions are out of date and run only the code to produce those targets.

The pipeline is invoked `targets::tar_make()`.  That runs the `_targets.R` file, which in turn sources `project_options.R` and the files within the `code` directory and, by default, produces all of the non-existent or out of date targets.

# Note

There may still be some aspects of this that have hardcoded variable names or
values, so if something doesn't work or looks funny I would first check for
those.

# Author

Miles Smith <miles-smith@omrf.org>

# Repository

https://github.com/milescsmith/rnaseq_targets_pipeline

## Last Updated
See the [CHANGELOG](CHANGELOG.md) for the latest updates
