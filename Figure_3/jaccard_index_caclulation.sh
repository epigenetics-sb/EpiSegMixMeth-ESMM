
# working directory
cd ./project/Analysis 

# conda core environment


# Not elegant way but using the following command can we calculate the jaccard index for all possible combinations for each state  
while read rfile; 
	do epirr=$(basename $rfile .bed.gz|cut -f 2 -d '_');
	ffile=$(grep "$epirr" ../files_list/EpiSegMix_meth_segmentation_bed_files.txt ) ;
	for i in $(seq 1 10);
		do for j in $(seq 1 10); 
			do si=$(bedtools jaccard -a <(zcat $ffile |awk -vOFS="\t" -v a="$i" '$4==a{print $0}' |sort -k1,1 -k2,2n -k3,3n) \
						 -b <(zcat $rfile|awk -vOFS="\t" -v b="$j" '$4==b{print $0}' |sort -k1,1 -k2,2n -k3,3n ) \
						 |tail -n+2 |cut -f 3); 		
						 echo $epirr $i $j $si \
						 |awk -vOFS="\t" '{print $1,$2,$3,sprintf("%f",$4)}' ;
			done ;
		done  ; 
	done < ../files_list/EpiSegMix_meth_segmentation_bed_files_4marks.txt |tee  jaccard_index_full_reduced/EpiSegMixMeth_full_reduced_si.txt
