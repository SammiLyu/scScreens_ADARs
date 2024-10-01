#!/bin/bash

mkdir -p collapsed_bam_H1_20230719
mkdir -p H1_parse_tmp

mkdir dedupped_bam_H1_20230721

mkdir -p AEI_H1_v2
mkdir -p AEI_H1_v2/tmp_20230126
mkdir -p AEI_H1_v2/log_20230126
mkdir -p AEI_H1_v2/index_20230126
mkdir -p AEI_H1_v2/mpileup_20230126

for FILE in cb_listsH1/v5/*
do

ter_id="$(echo $FILE | cut -d'_' -f3)" 
ter_id="$(echo $ter_id | cut -d'.' -f1)"
file_name="$(echo $FILE | cut -d'/' -f3)"
file_name="$(echo $file_name | cut -d'.' -f1)"
ter_id_num="$(echo $ter_id | cut -d'r' -f2)"


echo $FILE
echo $ter_id
echo $file_name
echo $ter_id_num

export BAM_FILE="/media/NAS1/Sammi_NAS1/CellRanger_Output/4H1/$ter_id_num/dm-$ter_id/outs/possorted_genome_bam.bam"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > H1_parse_tmp/SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > H1_parse_tmp/filtered_SAM_body
cat SAM_header H1_parse_tmp/filtered_SAM_body > H1_parse_tmp/filtered.sam
samtools view -b H1_parse_tmp/filtered.sam > collapsed_bam_H1_20230719/"$file_name".bam

samtools index collapsed_bam_H1_20230719/"$file_name".bam

done

for FILE in collapsed_bam_H1_20230719/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)_dedup.bam"

echo $file_name
echo $new_file_name

umi_tools dedup -I collapsed_bam_H1_20230719/$file_name -S dedupped_bam_H1_20230721/$new_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag
done


for FILE in dedupped_bam_H1_20230721/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)"
aei_input_file="$(echo $new_file_name).headerfixed.bam"

echo $file_name
echo $new_file_name
echo $aei_input_file


# only need to fix header for sc
samtools view -H dedupped_bam_H1_20230721/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - dedupped_bam_H1_20230721/$file_name > AEI_H1_v2/tmp_20230126/$new_file_name.headerfixed.bam

done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_H1_v2/tmp_20230126 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_H1_v2/log_20230126/ \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_H1_v2/output_index_20230126/ \
-f _dedup.headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_H1_v2/output_mpileup_20230126/ \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
-mm AllMismatches \
--verbose \
--stranded