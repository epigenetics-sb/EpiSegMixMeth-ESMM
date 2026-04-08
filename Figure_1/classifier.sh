#!/bin/bash
# Taking all the parameters....

while getopts ":hf:s:o:m:p:n:" OPTION
do
        case $OPTION in
                h)
                        printHelp; exit 0 ;;
                f)
                        segment_bed="$OPTARG" ;;
                s)
                        states=$OPTARG ;;
                o)
                        output_dir=$OPTARG ;;
		m)
                        mark=$OPTARG ;;
                p)
                        pilot_regions=$OPTARG ;;
                n) 
                        name=$OPTARG ;;
        esac
done

pname=${name}_${mark}_${states}N
temp_bed=${output_dir}/${pname}_temp.bed
segment_pilot=${output_dir}/${pname}_segment_pilot.bed
overlap_output=${output_dir}/${pname}_overlapEnrichment
ne_output=${output_dir}/${pname}_neighEnrichment
gene_output=${output_dir}/${pname}_geneExpression.txt
his_meth_overlap=${output_dir}/${pname}_heatmap_meth.txt
histone_overlap=${output_dir}/${pname}_heatmap.txt
heatmap_overlap=${output_dir}/${pname}_histone_meth.CpG.txt

# Function to process bed file based on mark
process_bed() {
    local mark=$1
    case $mark in
        EpiSegMix|EpiSegMixMeth)
            zcat ${segment_bed} | awk -v OFS="\t" '{print $0}' > ${temp_bed}
            ;;
        ChromHMM)
            cat ${segment_bed} | awk -v OFS="\t" 'NR>1{print $0}' > ${temp_bed}
            ;;
        EpiCSeg)
            cat ${segment_bed} | awk -v OFS="\t" '{print "chr"$0}' > ${temp_bed}
            ;;
        Segway)
            zcat ${segment_bed} | awk -v OFS="\t" 'NR>1{print $0}' > ${temp_bed}
            ;;
    esac
}

# Process the bed file
process_bed ${mark}

# Variables 
# coordinate files and anchor files for enrichment
anchor="./resources/references/ANCHORFILES/hg38/"

# meth tab file information and counts file 
methtab=$(grep ${name} ./project/files_list/wgbs_tab.txt)
countfile=./project/counts/EpiSegMix/${name}_refined_counts.txt
#countfile=./project/counts/EpiSegMix/${name}_refined_counts.v2.txt

# Filtering segment bedfile only for pilot regions
# segment_pilot=${temp_bed} # uncomment this line and comment the next if 
# bedtools intersect -a ${temp_bed} -b <(awk -vOFS="\t" '{print "chr"$0}' ${pilot_regions}) -u > ${segment_pilot} 

# Run the overlap enrichment for different coordinate sets
coord_paths=(
    "./resources/references/COORDS/hg38/"
    "./resources/references/LADs/bed/COORDS/"
    "./project/Analysis/cCREs/regions_hg38/"
    "./resources/references/temp/enhancers.bed"
)
output_suffixes=(
    "overlapEnrichment"
    "overlapEnrichment_LADs"
    "overlapEnrichment_cCREs"
    "overlapEnrichment_enhancers"
)

for i in "${!coord_paths[@]}"; do
    coords="${coord_paths[$i]}"
    overlap_output="${output_dir}/${pname}_${output_suffixes[$i]}"
    singularity exec -B ./mounted/ ./data/projects/container/Mixture_Model/chromhmm/ \
        java -Djava.awt.headless=true -jar -mx6000M -jar /ChromHMM/ChromHMM.jar OverlapEnrichment \
        -color 63,127,147 \
        ${segment_bed} \
        ${coords} \
        ${overlap_output}
done

#######################################################################
# !!!! COMMENT OR REMOVE THIS for NORMAL run
# exit 0

#######################################################################
<< 'COMMENT'
COMMENT
# Run the neighborhood enrichment
# using ChromHMM singularity containers 
singularity exec -B ./mounted/ ./data/projects/container/Mixture_Model/chromhmm/ \
    java -Djava.awt.headless=true -mx6000M -jar /ChromHMM/ChromHMM.jar  NeighborhoodEnrichment \
    -color 63,127,147  \
    ${temp_bed} \
    ${anchor}RefSeqTSS.hg38.txt.gz \
    ${ne_output}_TSS

singularity exec -B ./mounted/ ./data/projects/container/Mixture_Model/chromhmm/ \
    java -Djava.awt.headless=true -mx6000M -jar /ChromHMM/ChromHMM.jar  NeighborhoodEnrichment \
    -color 63,127,147  \
    ${temp_bed} \
    ${anchor}RefSeqTES.hg38.txt.gz \
    ${ne_output}_TES



touch ${gene_output}
for i in $(seq 1 ${states}); 
do
    gene_fpkm=$(cut -f 5,9,10 ${gene_overlap} |awk -v state="$i" '$2==state{print $1}' |sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }') ; \
    echo -e "$i\t$gene_fpkm" >> ${gene_output}
done

# DNA methylation, CpG density and histone marks overlap 
source ./env/miniconda3/bin/activate ./env/miniconda3/envs/core

nhistone=$(echo $header |wc -w)


header=$(head -n 1 ${countfile}|cut -f 4-)

echo "Generating the histone counts"
bedtools intersect -a ${temp_bed} -b <(awk -vOFS=\t 'NR>1{print "chr"$0}' ${countfile}) -wa -wb \
        |cut -f 1-4,13-18 \
        |bedtools groupby -i - -g 1,2,3,4 -c 5,6,7,8,9,10 -o median \
        > ${histone_overlap}

echo "Adding methylation counts"
bedtools intersect -a  ${histone_overlap} -b <(awk -vOFS="\t" '{print $1,$2,$2+1,$3,$4}' ${methtab}) -wa -wb \
        |bedtools groupby  -g 1,2,3,4 -c 5,6,7,8,9,10,14,15 -o first,first,first,first,first,first,mean,mean \
        |awk -v a="$header" -vOFS="\t" 'BEGIN{print "chr\tstart\tend\tstate\t"a"\tCov\tMeth"}{print $0}' \
        > ${his_meth_overlap} 

echo "Adding CpG counts"
bedtools intersect -a  <(awk -vOFS="\t" 'NR>1{print $0}' ${his_meth_overlap}) -b <(awk -vOFS="\t" '{print $1,$2,$2+1,$3,$4}' ${methtab}) -c \
        |sort -k1,1 -k2,2n -k3,3n \
        |awk -vOFS="\t" -v a="$header" 'BEGIN{print "chr\tstart\tend\tstate\t"a"\tCov\tMeth\tCpG"}{print $0}' \
        > ${heatmap_overlap}

echo "Generating the average histone, methylation and CpG counts for each state"
Rscript --vanilla ./project/src/classifier.r  ${heatmap_overlap} ${mark}_${name} ${output_dir}/${pname} ${states} ${nhistone}  &> ${output_dir}/${pname}_Rcommand.out

# cleaning dir
rm $histone_overlap $his_meth_overlap $temp_bed
echo "Classifier script finished!" 