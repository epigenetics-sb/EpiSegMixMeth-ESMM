#!/bin/bash
#SBATCH --job-name="Classifier_${SLURM_ARRAY_TASK_ID}"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-154%4
#SBATCH --mem=40G
#SBATCH --output=./project/slurm_logs/Classifier_%A_%a.out
#SBATCH --error=./project/slurm_logs/Classifier_%A_%a.err

# Specify the path to the config file
config="./project/counts/config.txt"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID and its variables
samplesheet=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
paired_end=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Define paths
pilot_regions="./data/episegmix/src/encode_pilot_regions/hg38.bed"
classifier="./project/src/classifier.sh"
states=10
source ./env/miniconda3/bin/activate ./env/miniconda3/envs/core

# Set mode: "normal" for reduced models, "4marks" for 4 marks
# MODE="4marks"
MODE="normal"

# Set bed file lists based on mode
if [ "$MODE" = "4marks" ]; then
    EpiSegMix_all_bed="./project/files_list/EpiSegMix_segmentation_bed_files_4marks.txt"
    EpiSegMix_meth_all_bed="./project/files_list/EpiSegMix_meth_segmentation_bed_files_4marks.txt"
    ChromHMM_all_bed="./project/files_list/ChromHMM_segmentation_bed_files_4marks.txt"
else
    EpiSegMix_all_bed="./data/shared/IHEC/Analysis/heatmap_matrix/EpiSegMix_segmentation_bed_files_reduced_models_H3K9me3.txt"
    EpiSegMix_meth_all_bed="./data/shared/IHEC/Analysis/heatmap_matrix/EpiSegMixMeth_segmentation_bed_files_reduced_models_H3K9me3.txt"
    ChromHMM_all_bed="./data/shared/IHEC/Analysis/heatmap_matrix/ChromHMM_segmentation_bed_files_reduced_models_H3K9me3.txt"
fi

# Set output base based on mode
if [ "$MODE" = "4marks" ]; then
    base_output="./project/results"
else
    base_output="./data/shared/IHEC/Analysis/heatmap_matrix/${model}"
fi

# Function to run classifier
run_classifier() {
    local bed_file=$1
    local output_dir=$2
    local model_name=$3
    bash ${classifier} -f ${bed_file} -s ${states} -o ${output_dir} -m ${model_name} -p ${pilot_regions} -n ${sample}
}

if [ "$MODE" != "4marks" ]; then
    mkdir -p ${base_output}
fi

# Process each model
for model_name in EpiSegMix EpiSegMixMeth ChromHMM; do
    if [ "$MODE" = "4marks" ]; then
        output_dir="${base_output}/${model_name}_4_marks/${states}N/${sample}"
        mkdir -p ${output_dir}
        if [ "$model_name" = "ChromHMM" ]; then
            bed_file="./project/results/ChromHMM_4marks/${states}N/${sample}/${sample}_${states}_dense.bed"
        else
            bed_list_var="${model_name}_all_bed"
            bed_list=${!bed_list_var}
            bed_file=$(grep ${sample} ${bed_list})
        fi
    else
        output_dir=${base_output}
        bed_list_var="${model_name}_all_bed"
        bed_list=${!bed_list_var}
        bed_file=$(grep ${sample} ${bed_list})
    fi
    echo "Processing ${model_name} for sample ${sample} in ${MODE} mode"
    run_classifier ${bed_file} ${output_dir} ${model_name}
done


# Later merge all results into one file for all classifers from all EpiRRs manually using bash to use the input for the Random Forest classifier resulting in the final files. 
# IHEC_ChromHMM_10N_all_states_samples_classifier.csv
# IHEC_EpiSegMix_10N_all_states_samples_classifier.csv
# IHEC_EpiSegMix_meth_10N_all_states_samples_classifier.csv

# One will also need manually labelled states for RFC. Use this to train to obtain the final labels across all samples and all classifiers. This will be the final output of the project.
# DEEP_samples_annotation_summary.tab
# DEEP_Bcells_trained_EpiRR_annotation_summary_chmm.csv
# DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix.csv
# DEEP_Bcells_trained_EpiRR_annotation_summary_episegmix_meth.csv


