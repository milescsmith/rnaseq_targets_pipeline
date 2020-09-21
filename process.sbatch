#! /bin/bash -l

#SBATCH -J sle_wb_rnaseq
#SBATCH -o all_sle_wb.log
#SBATCH --mail-user=miles-smith@omrf.org
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=256
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32

module load singularity
singularity run \
    --bind /s/guth-aci/all_sle_wb_rnaseq_samples:/home/rstudio/analysis/rnaseq/all_sle_wb_rnaseq_samples,/s/guth-aci/novaseq/:/home/rstudio/datasets/rnaseq/novaseq \
    /Volumes/guth_aci_informatics/software/containers/workerbee-ml-v1.1.sif \
    Rscript -e "drake::r_make('/home/rstudio/analysis/rnaseq/all_sle_wb_rnaseq_samples/_drake.R')"
