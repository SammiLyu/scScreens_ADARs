#!/bin/bash

input_dir=/media/NAS1/Sammi_NAS1/CellRanger_Output/cortical_organoids_cr

mkdir -p collapsed_bam_CO_20230627
mkdir -p collapsed_bam_CO_20230627/tmp
mkdir -p tagged_bam_CO_20230627
mkdir -p dedupped_bam_CO_20230627

mkdir AEI_CO_v2
mkdir AEI_CO_v2/tmp_20230627
mkdir AEI_CO_v2/log_20230627
mkdir AEI_CO_v2/output_index_20230627
mkdir AEI_CO_v2/output_mpileup_20230627

for FILE in cb_listsCO/v2/*
do
ter_id="$(echo $FILE | cut -d'/' -f3)" 
ter_id="$(echo $ter_id | cut -d'_' -f1)"
file_name="$(echo $FILE | cut -d'/' -f3)"
file_name="$(echo $file_name | cut -d'.' -f1)"
file_existence="/media/NAS1/Sammi_NAS1/CellRanger_Output/cortical_organoids_cr/$(echo $ter_id)/co_$(echo $ter_id)/outs/possorted_genome_bam.bam"
parsed_bam="collapsed_bam_CO_20230627/$(echo $file_name ).bam"

echo $FILE
echo $ter_id
echo $file_name
echo $file_existence

#if [ ! -f "$file_sorted_bam" ]; then
#samtools sort $file_existence -o $file_sorted_bam
#fi

if [ ! -f "$parsed_bam" ]; then
export BAM_FILE="$file_existence"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > collapsed_bam_CO_20230627/tmp/SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > collapsed_bam_CO_20230627/tmp/filtered_SAM_body
cat collapsed_bam_CO_20230627/tmp/SAM_header collapsed_bam_CO_20230627/tmp/filtered_SAM_body > collapsed_bam_CO_20230627/tmp/filtered.sam
samtools view -b collapsed_bam_CO_20230627/tmp/filtered.sam > $parsed_bam

samtools index $parsed_bam
fi

done

rm -r collapsed_bam_CO_20230627/tmp

for FILE in collapsed_bam_CO_20230627/*.bam

do
in_file_name="$(echo $FILE | cut -d'/' -f2)"
out_file_name="$(echo $in_file_name | cut -d'.' -f1)_dedup.bam"

echo $in_file_name
echo $out_file_name

if [ ! -f "dedupped_bam_CO_20230627/$dedupped_file_name" ]; then
umi_tools dedup -I collapsed_bam_CO_20230627/$in_file_name -S dedupped_bam_CO_20230627/$out_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag
fi
done

for FILE in dedupped_bam_CO_20230627/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)"

echo $file_name
echo $new_file_name

# only need to fix header for sc
samtools view -H dedupped_bam_CO_20230627/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - dedupped_bam_CO_20230627/$file_name > AEI_CO_v2/tmp_20230627/$new_file_name.headerfixed.bam

done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_CO_v2/tmp_20230627 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_CO_v2/log_20230627/AEI_CO_v2_log \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_CO_v2/output_index_20230627/AEI_CO_v2_outputindex \
-f _dedup.headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_CO_v2/output_mpileup_20230627/AEI_CO_v2_outputmpileup \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
-mm AllMismatches \
--stranded \
--verbose  
