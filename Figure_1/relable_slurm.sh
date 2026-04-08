#!/bin/bash
#SBATCH --job-name="Relable_${SLURM_ARRAY_TASK_ID}"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-154%4
#SBATCH --mem=40G
#SBATCH --output=./project/slurm_logs/Relabel_%A_%a.out
#SBATCH --error=./project/slurm_logs/Relabel_%A_%a.err

# Specify the path to the config file
config="./project/counts/config.txt"
#config="./project/IHEC/counts/config_ge_failed.txt"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID and its variables 
samplesheet=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

paired_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)


# set mode: episegmix | episegmix_meth | chmm
MODE="episegmix"

states=10
source ./env/miniconda3/bin/activate ./env/miniconda3/envs/core
# extract the chromHMM, EpiCSeg and EpiSegMix (hsitone mark model) as well as gene expression fpkm file 

relabel_script="./project/src/relable.r"
base_dir="./resources/IHEC"
base_out="${base_dir}/Analysis/relabelled_bed_files/random_forest_classifier"
base_ann="${base_dir}/Analysis/classifier/Random_forest_classifier"

case "${MODE}" in
  episegmix)
    bed=$(grep ${sample} "${base_dir}/files_list/EpiSegMix_segmentation_bed_files.txt")
    output="${base_out}/EpiSegMix/"
    annotation="${base_ann}/EpiSegMix_annotation_without_features.csv"
    tag="EpiSegMix"
    ;;
  episegmix_meth)
    bed=$(grep ${sample} "${base_dir}/files_list/EpiSegMix_meth_segmentation_bed_files.txt")
    output="${base_out}/EpiSegMixMeth/"
    annotation="${base_ann}/EpiSegMixMeth_annotation_without_features.csv"
    tag="EpiSegMixMeth"
    ;;
  chmm)
    bed="${base_dir}/results/ChromHMM_${states}N/${sample}/${sample}_${states}_dense.bed"
    output="${base_out}/ChromHMM/"
    annotation="${base_ann}/ChromHMM_annotation_without_features.csv"
    tag="ChromHMM"
    ;;
  *)
    echo "Unknown MODE: ${MODE}. Use episegmix | episegmix_meth | chmm" >&2
    exit 1
    ;;
esac

Rscript ${relabel_script} ${bed} ${output} ${annotation} ${sample} ${tag}