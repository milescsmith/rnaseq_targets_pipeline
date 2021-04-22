{targets}-based RNAseq analysis pipeline

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
