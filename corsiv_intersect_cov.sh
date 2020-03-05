for file in ./USC_pediatric_glioma_alignments/merged_CpG_eviddence/*.cov
do
	echo $(echo $(basename $file) | cut -f1 -d ".").bed
	bedtools intersect -a $file -b combined_corsivs.bed -wb > ./USC_pediatric_glioma_bed_file/$(echo $(basename $file) | cut -f1 -d ".").bed
done

