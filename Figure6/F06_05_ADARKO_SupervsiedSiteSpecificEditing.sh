#!/bin/bash

rm -r collapsed_bam_KO_20230801
mkdir collapsed_bam_KO_20230801
mkdir -p KO_tmp2

mkdir dedupped_bam_KO_20230801

mkdir edit_output_KO_v7

for FILE in cb_listsKO/v7/*.txt
do

ter_id="$(echo $FILE | cut -d'_' -f2)" 
ter_id="$(echo $ter_id | cut -d'/' -f3)" 
file_name="$(echo $FILE | cut -d'/' -f3)"
file_name="$(echo $file_name | cut -d'.' -f1)"
#ter_id_num="$(echo $ter_id | cut -d'r' -f2)"

echo $FILE
echo $ter_id
echo $file_name
#echo $ter_id_num

export BAM_FILE="/media/NAS1/Sammi_NAS1/Cellranger_Output/ADAR_20220614/dual_ref/$ter_id/outs/possorted_genome_bam.bam"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > KO_tmp2/SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > KO_tmp2/filtered_SAM_body
cat KO_tmp2/SAM_header KO_tmp2/filtered_SAM_body > KO_tmp2/filtered.sam
samtools view -b KO_tmp2/filtered.sam > collapsed_bam_KO_20230801/"$file_name".bam

samtools index collapsed_bam_KO_20230801/"$file_name".bam

done

for FILE in collapsed_bam_KO_20230801/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)_dedup.bam"

echo $file_name
echo $new_file_name

umi_tools dedup -I collapsed_bam_KO_20230801/$file_name -S dedupped_bam_KO_20230801/$new_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag

done

for FILE in dedupped_bam_KO_20230801/*.bam
do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1).txt"

echo $file_name
echo $new_file_name

perl query_known_sites.pl /media/Scratch_SSD_Voyager/sammi/RNA_editing/TABLE1_GRCh38.txt /media/Scratch_SSD_Voyager/sammi/RNA_editing/dedupped_bam_KO_20230801/$file_name edit_output_KO_v7/$new_file_name
done