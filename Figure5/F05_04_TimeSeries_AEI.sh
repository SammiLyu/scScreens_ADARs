#!/bin/bash

rm -r collapsed_bam_TS_20230124
mkdir collapsed_bam_TS_20230124
mkdir dedupped_bam_TS_20230125

mkdir -p AEI_TS_v4rerunverbose
mkdir -p AEI_TS_v4rerunverbose/tmp_20230623
mkdir -p AEI_TS_v4rerunverbose/log_20230623
mkdir -p AEI_TS_v4rerunverbose/output_index_20230623
mkdir -p AEI_TS_v4rerunverbose/output_mpileup_20230623


for FILE in cb_listsTS/v1/*
do

ter_id="$(echo $FILE | cut -d'_' -f3)_$(echo $FILE | cut -d'_' -f4)_$(echo $FILE | cut -d'_' -f5)"
ter_id="$(echo $ter_id | cut -d'.' -f1)"
file_name="$(echo $FILE | cut -d'/' -f3)"
file_name="$(echo $file_name | cut -d'.' -f1)"
ter_id_num="$(echo $ter_id | cut -d'r' -f2)"


echo $FILE
echo $ter_id
echo $file_name

RNA="wk4_M1_1 wk6_M15_1 wk8_M8_1 wk10_A_2 wk4_M14_1 wk6_M9_1 wk8_M10_1 wk10_B_1"

[[ $RNA =~ (^|[[:space:]])$ter_id($|[[:space:]]) ]] && export BAM_FILE="/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/Teratoma_TS_RNA_intronT/multispecies/$ter_id/outs/possorted_genome_bam.bam" || export BAM_FILE="/media/NAS1/Sammi_NAS1/CellRanger_Output/MLT_20220103_TS_NovaSeq_10xARC/$ter_id/outs/gex_possorted_bam.bam"

echo $BAM_FILE

source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > filtered_SAM_body
cat SAM_header filtered_SAM_body > filtered.sam
samtools view -b filtered.sam > collapsed_bam_TS_20230124/"$file_name".bam

samtools index collapsed_bam_TS_20230124/"$file_name".bam

done

for FILE in collapsed_bam_TS_20230124/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)_dedup.bam"

echo $file_name
echo $new_file_name

umi_tools dedup -I collapsed_bam_TS_20230124/$file_name -S dedupped_bam_TS_20230125/$new_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag

done

for FILE in dedupped_bam_TS_20230125/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)"

echo $file_name
echo $new_file_name

# only need to fix header for sc
samtools view -H dedupped_bam_TS_20230125/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - dedupped_bam_TS_20230125/$file_name > AEI_TS_v4rerunverbose/tmp_20230623/$new_file_name.headerfixed.bam

done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_TS_v4rerunverbose/tmp_20230623 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_TS_v4rerunverbose/log_20230623/AEI_TS_v4rerunverbose_log \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_TS_v4rerunverbose/output_index_20230623/AEI_TS_v4rerunverbose_outputindex \
-f _dedup.headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_TS_v4rerunverbose/output_mpileup_20230623/AEI_TS_v4rerunverbose_outputmpileup \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
-mm AllMismatches \
--stranded \
--verbose

