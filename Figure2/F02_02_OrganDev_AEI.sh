#!/bin/bash

# Data downloaded from https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6814

input_dir=/media/NAS1/Sammi_NAS1/Processed_Reads/OrganDevelopment/human
samtools_tmp_dir=OrganDev_v1

mkdir -p $samtools_tmp_dir/fixmate
mkdir -p $samtools_tmp_dir/possorted
mkdir -p $samtools_tmp_dir/markdup
mkdir -p $samtools_tmp_dir/nomultiseg

mkdir -p AEI_OrganDev_v2
mkdir -p AEI_OrganDev_v2/tmp_20230607
mkdir -p AEI_OrganDev_v2/log_20230607
mkdir -p AEI_OrganDev_v2/output_index_20230607
mkdir -p AEI_OrganDev_v2/output_mpileup_20230607

for FILE in $input_dir/*.bam
do

file_name="$(echo $FILE | cut -d'/' -f8)"
sample="$(echo $file_name | cut -d'.' -f1).$(echo $file_name | cut -d'.' -f2).$(echo $file_name | cut -d'.' -f3).$(echo $file_name | cut -d'.' -f4).$(echo $file_name | cut -d'.' -f5)"
fixmate_file_name="$(echo $sample).fixmate.bam"
possorted_file_name="$(echo $sample).possorted.bam"
markdup_file_name="$(echo $sample).markdup.bam"
headerfix_file_name="$(echo $sample).headerfixed.bam"
nomultiseg_file_name="$(echo $sample).nomultiseg.bam"

# bam is already name sorted so skip name sorting
#samtools sort -n -o $samtools_tmp_dir/namesort/$namesort_file_name $samtools_tmp_dir/$file_name

if [ ! -f "$samtools_tmp_dir/fixmate/$fixmate_file_name" ]; then
    echo "running samtools fixmate on "$sample
    samtools fixmate -m $input_dir/$file_name $samtools_tmp_dir/fixmate/$fixmate_file_name -@ 4
fi

if [ ! -f "$samtools_tmp_dir/possorted/$possorted_file_name" ]; then
    echo "running samtools position sort on "$sample
    samtools sort -o $samtools_tmp_dir/possorted/$possorted_file_name $samtools_tmp_dir/fixmate/$fixmate_file_name -@ 4
fi

if [ ! -f "$samtools_tmp_dir/markdup/$markdup_file_name" ]; then
    echo "running remove duplicates on "$sample
    samtools markdup -r $samtools_tmp_dir/possorted/$possorted_file_name $samtools_tmp_dir/markdup/$markdup_file_name -@ 4
    printf "\n"
fi

# no need to fix header for this dataset at all
#if [ ! -f "AEI_$samtools_tmp_dir/tmp_20230607/$headerfix_file_name" ]; then
#    echo "fixing header for "$sample
    #samtools view -H $samtools_tmp_dir/markdup/$markdup_file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - $samtools_tmp_dir/markdup/$markdup_file_name > AEI_OrganDev_v2/tmp_20230607/$headerfix_file_name
#fi

# remove all multi segmented reads - avoid pair-end error
if [ ! -f "$samtools_tmp_dir/nomultiseg/$nomultiseg_file_name" ]; then
    echo "running remove multi-segmented reads on "$sample
    samtools view -F 1 -b $samtools_tmp_dir/markdup/$markdup_file_name > $samtools_tmp_dir/nomultiseg/$nomultiseg_file_name
    printf "\n"
fi

done

# check if all files in directory have been processed
ls $input_dir/*.bam | sed 's#/media/NAS1/Sammi_NAS1/Processed_Reads/OrganDevelopment/human#OrganDev_v1/nomultiseg#' | sed 's#sorted.bam#nomultiseg.bam#' > AEI_OrganDev_v2/aei_input_bam_summary_compto.txt
ls $samtools_tmp_dir/nomultiseg/*.bam > AEI_OrganDev_v2/aei_input_bam_summary.txt

if cmp -s "AEI_OrganDev_v2/aei_input_bam_summary_compto.txt" "AEI_OrganDev_v2/aei_input_bam_summary.txt"; then
    printf 'Finished processing all input bams. Moving on to AEI...\n'
    /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/$samtools_tmp_dir/nomultiseg \
    -l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_OrganDev_v2/log_20230607/OrganDev_v1_log \
    -os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_OrganDev_v2/output_index_20230607/OrganDev_v1_outputindex \
    -f .nomultiseg.bam  \
    -o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_OrganDev_v2/output_mpileup_20230607/OrganDev_v1_outputmpileup \
    --genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg19GTExGeneExpression.bed.gz \
    --refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg19RefSeqCurated.bed.gz \
    --snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg19CommonGenomicSNPs150.bed.gz \
    -gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg19Genome.fa \
    -rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg19Alu.bed.gz \
    --genome hg19 \
    -mm AllMismatches \
    --verbose \
    --stranded
else
    printf 'Not all bams are processed. User attention required.\n'
fi

if [ -f "AEI_OrganDev_v2/output_index_20230607/OrganDev_v1_outputindex/EditingIndex.csv" ]; then
    echo "Editing index table successfully generated. Removing intermediate processing files...\n"
    #rm -r $samtools_tmp_dir
fi

