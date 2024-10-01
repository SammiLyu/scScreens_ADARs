#!/bin/bash
rm -r collapsed_bam_H1_20230102
mkdir collapsed_bam_H1_20230102
mkdir dedupped_bam_H1_20230102

mkdir AEI_H1_v1
mkdir AEI_H1_v1/tmp_20230518
mkdir AEI_H1_v1/log_20230518
mkdir AEI_H1_v1/output_index_20230518
mkdir AEI_H1_v1/output_mpileup_20230518


for FILE in cb_listsH1/v4/*
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

export BAM_FILE="/media/Scratch_SSD_Voyager/sammi/4H1/$ter_id_num/dm-$ter_id/outs/possorted_genome_bam.bam"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > filtered_SAM_body
cat SAM_header filtered_SAM_body > filtered.sam
samtools view -b filtered.sam > collapsed_bam_H1_20230102/"$file_name".bam

samtools index collapsed_bam_H1_20230102/"$file_name".bam

done

for FILE in collapsed_bam_H1_20230102/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)_dedup.bam"

echo $file_name
echo $new_file_name

umi_tools dedup -I collapsed_bam_H1_20230102/$file_name -S dedupped_bam_H1_20230102/$new_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag
done

for FILE in dedupped_bam_H1_20230102/*.bam
do
file_name="$(echo $FILE | cut -d'/' -f2)"
new_file_name="$(echo $file_name | cut -d'.' -f1)"
aei_input_file="$(echo $new_file_name).headerfixed.bam"

echo $file_name
echo $new_file_name
echo $aei_input_file

# only need to fix header for sc
samtools view -H dedupped_bam_H1_20230102/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - dedupped_bam_H1_20230102/$file_name > AEI_H1_v1/tmp_20230518/$new_file_name.headerfixed.bam
done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_H1_v1/tmp_20230518 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_H1_v1/log_20230518/ \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_H1_v1/output_index_20230518/ \
-f _dedup.headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_H1_v1/output_mpileup_20230518/ \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
#-mm AllMismatches \
--stranded
