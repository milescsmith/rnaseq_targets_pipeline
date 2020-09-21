# CLE RNAseq data - 2020/09/08

## Goals
- Differential gene analysis of CLE (and subtypes) vs Control subjects
- Score previously described gene modules (Banchereau, et al. 2016 10.1016/j.cell.2016.03.008; Kegerreis, et al. 2019 10.4049/jimmunol.1801512; Haynes, et al. 2019  10.1101/834093)
- Identify patient clusters
- Identify gene expression modules associated with particular clusters, disease-classifications

## Analyst
- Miles Smith <miles-smith@omrf.org>

## Description
RNA-sequencing of RNA isolated from total PBMCs.  RNA isolated by Qiagen RNeasy, libraries contructed using
KAPA Hyperprep (with ribosomal and globin depletion), and sequenced on a Novaseq 6000 S4.  Raw data was
converted from bcls to FASTQs using bcl2fastq 2.20, trimmed using bbmap 38.72, pseudoaligned and quantified
with Salmon 1.0.0 using a index built from GENCODE v32.  All but the bcl -> FASTQ conversion was performed
using a Snakemake script (https://gitlab.com/milothepsychic/rnaseq_pipeline/,
commit df16b508bc118350c82ebc2dcd485d7cc9c59a5f) and docker containers.

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

## Data sources
RNA was extracted by Laura Battiest in the Biorepository and sequencing libraries
were created by Nicolas Dominguez and Tate Bolin.  Sequencing was perfomed in the
OMRF Clinical Genomics Core.  RNAseq data was processed by me (Miles Smith) as
detailed above.

Cytokine levels were generated on a Bioplex 200 from serum|plasma (?) by Tim Gross.

Clinical, ANA, and cytokine data was compiled from the file "Cutaneous Lupus Project_Dr.Chong_190805.xlsx"
found in Ly Tran's folders.

## Usage

Modify the `_drake.r` file so that variables in the `### Setup file locations` point to the correct locations.
In the directory containing this project, run:
  ```
docker run -it --rm --mount type=bind,source=<location_of_files>,target=<location_entered_in_drake.r> gitlab.com/guthridge_informatics/control/workerbee-rnaseq-drake:1.4 /bin/Rscript "drake::r_make()"
```

(Note: this analysis may require a lot of memory.  It has been sucessfully run with 52G, but I am unclear how much is actually required.)

### Analysis on the OMRF HPC

Since Docker is incompatible with use on a shared compute cluster (Docker requires elevated permissions that no sane admin would grant), it is necessary to use [Singularity](https://sylabs.io/guides/3.5/user-guide/).  To run:
  ```
module load singularity
singularity exec --bind=<location_entered_in_drake.r> /Volumes/guth_aci_informatics/software/workerbee-sc.sif Rscript -e "drake::r_make('<path to _drake.r>')"
```

## Last Updated
See the [CHANGELOG](CHANGELOG.md) for the latest updates
