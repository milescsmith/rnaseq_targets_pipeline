{targets}-based RNAseq analysis pipeline

# Warning

Some aspects of this still have hardcoded variable names that are a hold over from my previous {drake}-based[rnaseq_drake](https://github.com/milescsmith/rnaseq_drake) pipeline.  Consider this still something of a work-in-progress.

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

The order of the columns does not matter.  Extra columns also do not make a difference.
See the template file in the `metadata` subdirectory.
