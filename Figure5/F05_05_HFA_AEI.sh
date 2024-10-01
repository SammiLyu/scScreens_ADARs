#!/bin/bash
mkdir -p HFA_v3/trimmed_fastq_HFA_20230711
mkdir -p HFA_v3/aligned_bam_HFA_20230711
mkdir -p HFA_v3/alignedfiltered_bam_HFA_20230711
mkdir -p HFA_v3/collapsed_bam_HFA_20230711
mkdir -p HFA_v3/tagged_bam_HFA_20230711
mkdir -p HFA_v3/dedupped_bam_HFA_20230711
mkdir -p HFA_v3/tmp_20230711

mkdir -p AEI_HFA_v5/tmp_20230711
mkdir -p AEI_HFA_v5/log_20230711
mkdir -p AEI_HFA_v5/output_index_20230711
mkdir -p AEI_HFA_v5/output_mpileup_20230711

#STAR --genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ --genomeLoad Remove
for FILE in cb_listsHFA/v3/*
do
ter_id="$(echo $FILE | cut -d'_' -f4)" 
file_name="$(echo $FILE | cut -d'/' -f3)"
ter_id_save="$(echo $file_name | cut -d'_' -f1).$(echo $file_name | cut -d'_' -f2)"
cell_type="$(echo $file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
file_name="$(echo $file_name | cut -d'_' -f1)_$(echo $file_name | cut -d'_' -f2)_$(echo $file_name | cut -d'_' -f3)_$cell_type"
ter_id_num="$(echo $ter_id_save).fastq.gz"
trimmed_fastq="HFA_v3/trimmed_fastq_HFA_20230711/$(echo $ter_id_save)_trimmed.fq.gz"
aligned_bam="HFA_v3/aligned_bam_HFA_20230711/$(echo $ter_id_save)Aligned.sortedByCoord.out.bam"
alignedfiltered_bam="HFA_v3/alignedfiltered_bam_HFA_20230711/$(echo $ter_id_save)AlignedFiltered.sortedByCoord.out.bam"
file_existence="HFA_v3/alignedfiltered_bam_HFA_20230711/$(echo $ter_id_save)AlignedFiltered.sortedByCoord.out.bam"
parsed_bam="HFA_v3/collapsed_bam_HFA_20230711/$(echo $file_name ).bam"

#file_bam="/media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$(echo $ter_id_save)/$(echo $ter_id_save)Aligned.out.bam"
#file_sorted_bam="/media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$(echo $ter_id_save)/$(echo $ter_id_save)Aligned.possorted.out.bam"

#echo $FILE
#echo $ter_id
#echo $ter_id_num
echo $file_name
#echo $ter_id_save
#echo $file_existence

if [ ! -f "$trimmed_fastq" ]; then
echo "$trimmed_fastq does not exist."
echo "trimming fastq..."
/media/Home_Raid1_Voyager/sammi/trim_galore/TrimGalore-0.6.10/trim_galore /media/NAS1/Sammi_NAS1/HFA_Raw/$ter_id/$ter_id_num \
-a AAAAAAAA --three_prime_clip_R1 1 \
-o HFA_v3/trimmed_fastq_HFA_20230711
fi

#mkdir -p /media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$ter_id_save
if [ ! -f "$aligned_bam" ]; then
echo "$aligned_bam does not exist."
echo "aligning using STAR..."
STAR --runThreadN 6 --outSAMstrandField intronMotif \
--genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ \
--readFilesCommand zcat \
--readFilesIn $trimmed_fastq \
--outFileNamePrefix HFA_v3/aligned_bam_HFA_20230711/$ter_id_save \
--outStd SAM \
--outSAMmultNmax -1 \
--outSAMtype BAM SortedByCoordinate \
--limitOutSJcollapsed 5000000 \
--limitBAMsortRAM 25000000000
# --genomeLoad LoadAndKeep
fi

#if [ ! -f "$file_sorted_bam" ]; then
#samtools sort $file_existence -o $file_sorted_bam
#fi

if [ ! -f "$alignedfiltered_bam" ]; then
samtools view -bh -q 30 -F 4 $aligned_bam | samtools sort -@ 10 - | samtools view -bh -> $alignedfiltered_bam
fi

if [ ! -f "$parsed_bam" ]; then
export BAM_FILE="$alignedfiltered_bam"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > HFA_v3/tmp_20230711/SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > HFA_v3/tmp_20230711/filtered_SAM_body
cat HFA_v3/tmp_20230711/SAM_header HFA_v3/tmp_20230711/filtered_SAM_body > HFA_v3/tmp_20230711/filtered.sam
samtools view -b HFA_v3/tmp_20230711/filtered.sam > $parsed_bam

samtools index $parsed_bam
fi


done
#STAR --genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ --genomeLoad Remove

for FILE in HFA_v3/collapsed_bam_HFA_20230711/*.bam

do
in_file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $in_file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
out_file_name="$(echo $in_file_name | cut -d'_' -f1)_$(echo $in_file_name | cut -d'_' -f2)_$(echo $in_file_name | cut -d'_' -f3)_$cell_type"
out_file_name="$(echo $out_file_name)_tag.bam"

echo $in_file_name
echo $out_file_name

if [ ! -f "HFA_v3/tagged_bam_HFA_20230711/$out_file_name" ]; then
python HFA_tag_20230601.py -i HFA_v3/collapsed_bam_HFA_20230711/$in_file_name -o HFA_v3/tagged_bam_HFA_20230711/$out_file_name
fi

done

for FILE in HFA_v3/tagged_bam_HFA_20230711/*.bam
do

tagged_file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $tagged_file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
dedupped_file_name="$(echo $tagged_file_name | cut -d'_' -f1)_$(echo $tagged_file_name | cut -d'_' -f2)_$(echo $tagged_file_name | cut -d'_' -f3)_$cell_type"
dedupped_file_name="$(echo $dedupped_file_name)_dedup.bam"

echo $tagged_file_name
echo $dedupped_file_name


#umi_tools dedup -I collapsed_bam_HFA_20230711/H27472_Eye_SRR12355470_Astrocytes.bam -S dedupped_bam_HFA_20230711/H27472_Eye_SRR12355470_Astrocytes_dedup.bam --umi-separator=",([^;]*)" --method=unique --extract-umi-method=read_id

if [ ! -f "HFA_v3/dedupped_bam_HFA_20230711/$dedupped_file_name" ]; then
samtools index HFA_v3/tagged_bam_HFA_20230711/$tagged_file_name
/media/Home_Raid1_Voyager/sammi/.local/bin/umi_tools dedup -I HFA_v3/tagged_bam_HFA_20230711/$tagged_file_name -S HFA_v3/dedupped_bam_HFA_20230711/$dedupped_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag
fi
done

for FILE in HFA_v3/dedupped_bam_HFA_20230711/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
new_file_name="$(echo $file_name | cut -d'_' -f1)_$(echo $file_name | cut -d'_' -f2)_$(echo $file_name | cut -d'_' -f3)_$cell_type"

aei_input_file="$(echo $new_file_name).headerfixed.bam"

echo $file_name
echo $new_file_name
echo $aei_input_file

# only need to fix header for sc
samtools view -H HFA_v3/dedupped_bam_HFA_20230711/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - HFA_v3/dedupped_bam_HFA_20230711/$file_name > AEI_HFA_v5/tmp_20230711/$new_file_name.headerfixed.bam

done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v5/tmp_20230711 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v5/log_20230711/AEI_HFA_v5_log \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v5/output_index_20230711/AEI_HFA_v5_outputindex \
-f .headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_HFA_v5/output_mpileup_20230711/AEI_HFA_v5_outputmpileup \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
-mm AllMismatches \
--stranded \
--verbose


mkdir -p HFA_v4/trimmed_fastq_HFA_20230711
mkdir -p HFA_v4/aligned_bam_HFA_20230711
mkdir -p HFA_v4/alignedfiltered_bam_HFA_20230711
mkdir -p HFA_v4/collapsed_bam_HFA_20230711
mkdir -p HFA_v4/tagged_bam_HFA_20230711
mkdir -p HFA_v4/dedupped_bam_HFA_20230711
mkdir -p HFA_v4/tmp_20230711

mkdir -p AEI_HFA_v6/tmp_20230711
mkdir -p AEI_HFA_v6/log_20230711
mkdir -p AEI_HFA_v6/output_index_20230711
mkdir -p AEI_HFA_v6/output_mpileup_20230711

#STAR --genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ --genomeLoad Remove
for FILE in cb_listsHFA/v4/*
do
ter_id="$(echo $FILE | cut -d'_' -f4)" 
file_name="$(echo $FILE | cut -d'/' -f3)"
ter_id_save="$(echo $file_name | cut -d'_' -f1).$(echo $file_name | cut -d'_' -f2)"
cell_type="$(echo $file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
file_name="$(echo $file_name | cut -d'_' -f1)_$(echo $file_name | cut -d'_' -f2)_$(echo $file_name | cut -d'_' -f3)_$cell_type"
ter_id_num="$(echo $ter_id_save).fastq.gz"
trimmed_fastq="HFA_v4/trimmed_fastq_HFA_20230711/$(echo $ter_id_save)_trimmed.fq.gz"
aligned_bam="HFA_v4/aligned_bam_HFA_20230711/$(echo $ter_id_save)Aligned.sortedByCoord.out.bam"
alignedfiltered_bam="HFA_v4/alignedfiltered_bam_HFA_20230711/$(echo $ter_id_save)AlignedFiltered.sortedByCoord.out.bam"
file_existence="HFA_v4/alignedfiltered_bam_HFA_20230711/$(echo $ter_id_save)AlignedFiltered.sortedByCoord.out.bam"
parsed_bam="HFA_v4/collapsed_bam_HFA_20230711/$(echo $file_name ).bam"

#file_bam="/media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$(echo $ter_id_save)/$(echo $ter_id_save)Aligned.out.bam"
#file_sorted_bam="/media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$(echo $ter_id_save)/$(echo $ter_id_save)Aligned.possorted.out.bam"

#echo $FILE
#echo $ter_id
#echo $ter_id_num
echo $file_name
#echo $ter_id_save
#echo $file_existence

if [ ! -f "$trimmed_fastq" ]; then
echo "$trimmed_fastq does not exist."
echo "trimming fastq..."
/media/Home_Raid1_Voyager/sammi/trim_galore/TrimGalore-0.6.10/trim_galore /media/NAS1/Sammi_NAS1/HFA_Raw/$ter_id/$ter_id_num \
-a AAAAAAAA --three_prime_clip_R1 1 \
-o HFA_v4/trimmed_fastq_HFA_20230711
fi

#mkdir -p /media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$ter_id_save
if [ ! -f "$aligned_bam" ]; then
echo "$aligned_bam does not exist."
echo "aligning using STAR..."
STAR --runThreadN 6 --outSAMstrandField intronMotif \
--genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ \
--readFilesCommand zcat \
--readFilesIn $trimmed_fastq \
--outFileNamePrefix HFA_v4/aligned_bam_HFA_20230711/$ter_id_save \
--outStd SAM \
--outSAMmultNmax -1 \
--outSAMtype BAM SortedByCoordinate \
--limitOutSJcollapsed 5000000 \
--limitBAMsortRAM 25000000000
# --genomeLoad LoadAndKeep
fi

#if [ ! -f "$file_sorted_bam" ]; then
#samtools sort $file_existence -o $file_sorted_bam
#fi

if [ ! -f "$alignedfiltered_bam" ]; then
samtools view -bh -q 30 -F 4 $aligned_bam | samtools sort -@ 10 - | samtools view -bh -> $alignedfiltered_bam
fi

if [ ! -f "$parsed_bam" ]; then
export BAM_FILE="$alignedfiltered_bam"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > HFA_v4/tmp_20230711/SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > HFA_v4/tmp_20230711/filtered_SAM_body
cat HFA_v4/tmp_20230711/SAM_header HFA_v4/tmp_20230711/filtered_SAM_body > HFA_v4/tmp_20230711/filtered.sam
samtools view -b HFA_v4/tmp_20230711/filtered.sam > $parsed_bam

samtools index $parsed_bam
fi


done
#STAR --genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ --genomeLoad Remove

for FILE in HFA_v4/collapsed_bam_HFA_20230711/*.bam

do
in_file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $in_file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
out_file_name="$(echo $in_file_name | cut -d'_' -f1)_$(echo $in_file_name | cut -d'_' -f2)_$(echo $in_file_name | cut -d'_' -f3)_$cell_type"
out_file_name="$(echo $out_file_name)_tag.bam"

echo $in_file_name
echo $out_file_name

if [ ! -f "HFA_v4/tagged_bam_HFA_20230711/$out_file_name" ]; then
python HFA_tag_20230601.py -i HFA_v4/collapsed_bam_HFA_20230711/$in_file_name -o HFA_v4/tagged_bam_HFA_20230711/$out_file_name
fi

done

for FILE in HFA_v4/tagged_bam_HFA_20230711/*.bam
do

tagged_file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $tagged_file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
dedupped_file_name="$(echo $tagged_file_name | cut -d'_' -f1)_$(echo $tagged_file_name | cut -d'_' -f2)_$(echo $tagged_file_name | cut -d'_' -f3)_$cell_type"
dedupped_file_name="$(echo $dedupped_file_name)_dedup.bam"

echo $tagged_file_name
echo $dedupped_file_name


#umi_tools dedup -I collapsed_bam_HFA_20230711/H27472_Eye_SRR12355470_Astrocytes.bam -S dedupped_bam_HFA_20230711/H27472_Eye_SRR12355470_Astrocytes_dedup.bam --umi-separator=",([^;]*)" --method=unique --extract-umi-method=read_id

if [ ! -f "HFA_v4/dedupped_bam_HFA_20230711/$dedupped_file_name" ]; then
samtools index HFA_v4/tagged_bam_HFA_20230711/$tagged_file_name
/media/Home_Raid1_Voyager/sammi/.local/bin/umi_tools dedup -I HFA_v4/tagged_bam_HFA_20230711/$tagged_file_name -S HFA_v4/dedupped_bam_HFA_20230711/$dedupped_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag
fi
done

for FILE in HFA_v4/dedupped_bam_HFA_20230711/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
new_file_name="$(echo $file_name | cut -d'_' -f1)_$(echo $file_name | cut -d'_' -f2)_$(echo $file_name | cut -d'_' -f3)_$cell_type"

aei_input_file="$(echo $new_file_name).headerfixed.bam"

echo $file_name
echo $new_file_name
echo $aei_input_file

# only need to fix header for sc
samtools view -H HFA_v4/dedupped_bam_HFA_20230711/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - HFA_v4/dedupped_bam_HFA_20230711/$file_name > AEI_HFA_v6/tmp_20230711/$new_file_name.headerfixed.bam

done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v6/tmp_20230711 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v6/log_20230711/AEI_HFA_v6_log \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v6/output_index_20230711/AEI_HFA_v6_outputindex \
-f .headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_HFA_v6/output_mpileup_20230711/AEI_HFA_v6_outputmpileup \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
-mm AllMismatches \
--stranded \
--verbose

mkdir -p HFA_v5/trimmed_fastq_HFA_20230711
mkdir -p HFA_v5/aligned_bam_HFA_20230711
mkdir -p HFA_v5/alignedfiltered_bam_HFA_20230711
mkdir -p HFA_v5/collapsed_bam_HFA_20230711
mkdir -p HFA_v5/tagged_bam_HFA_20230711
mkdir -p HFA_v5/dedupped_bam_HFA_20230711
mkdir -p HFA_v5/tmp_20230711

mkdir -p AEI_HFA_v7/tmp_20230711
mkdir -p AEI_HFA_v7/log_20230711
mkdir -p AEI_HFA_v7/output_index_20230711
mkdir -p AEI_HFA_v7/output_mpileup_20230711

#STAR --genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ --genomeLoad Remove
for FILE in cb_listsHFA/v5/*
do
ter_id="$(echo $FILE | cut -d'_' -f4)" 
file_name="$(echo $FILE | cut -d'/' -f3)"
ter_id_save="$(echo $file_name | cut -d'_' -f1).$(echo $file_name | cut -d'_' -f2)"
cell_type="$(echo $file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
file_name="$(echo $file_name | cut -d'_' -f1)_$(echo $file_name | cut -d'_' -f2)_$(echo $file_name | cut -d'_' -f3)_$cell_type"
ter_id_num="$(echo $ter_id_save).fastq.gz"
trimmed_fastq="/media/NAS1/Sammi_NAS1/RNA_editing_Archive/HFA_v2_Archive/trimmed_fastq_HFA_20230711/$(echo $ter_id_save)_trimmed.fq.gz"
aligned_bam="/media/NAS1/Sammi_NAS1/RNA_editing_Archive/HFA_v2_Archive/aligned_bam_HFA_20230711/$(echo $ter_id_save)Aligned.sortedByCoord.out.bam"
alignedfiltered_bam="HFA_v5/alignedfiltered_bam_HFA_20230711/$(echo $ter_id_save)AlignedFiltered.sortedByCoord.out.bam"
file_existence="/media/NAS1/Sammi_NAS1/RNA_editing_Archive/HFA_v2_Archive/alignedfiltered_bam_HFA_20230711/$(echo $ter_id_save)AlignedFiltered.sortedByCoord.out.bam"
parsed_bam="HFA_v5/collapsed_bam_HFA_20230711/$(echo $file_name ).bam"

#file_bam="/media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$(echo $ter_id_save)/$(echo $ter_id_save)Aligned.out.bam"
#file_sorted_bam="/media/Home_Raid1_Voyager/sammi/sammi_Scratch/HFA/$(echo $ter_id_save)/$(echo $ter_id_save)Aligned.possorted.out.bam"

#echo $FILE
#echo $ter_id
#echo $ter_id_num
echo $file_name
#echo $ter_id_save
#echo $file_existence

#if [ ! -f "$trimmed_fastq" ]; then
#echo "$trimmed_fastq does not exist."
#echo "trimming fastq..."
#/media/Home_Raid1_Voyager/sammi/trim_galore/TrimGalore-0.6.10/trim_galore /media/NAS1/Sammi_NAS1/HFA_Raw/$ter_id/$ter_id_num \
#-a AAAAAAAA --three_prime_clip_R1 1 \
#-o HFA_v5/trimmed_fastq_HFA_20230711
#fi

#if [ ! -f "$aligned_bam" ]; then
#echo "$aligned_bam does not exist."
#echo "aligning using STAR..."
#STAR --runThreadN 6 --outSAMstrandField intronMotif \
#--genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ \
#--readFilesCommand zcat \
#--readFilesIn $trimmed_fastq \
#--outFileNamePrefix HFA_v5/aligned_bam_HFA_20230711/$ter_id_save \
#--outStd SAM \
#--outSAMmultNmax -1 \
#--outSAMtype BAM SortedByCoordinate \
#--limitOutSJcollapsed 5000000 \
#--limitBAMsortRAM 25000000000
# --genomeLoad LoadAndKeep
#fi


#if [ ! -f "$alignedfiltered_bam" ]; then
#samtools view -bh -q 30 -F 4 $aligned_bam | samtools sort -@ 10 - | samtools view -bh -> $alignedfiltered_bam
#fi

if [ ! -f "$parsed_bam" ]; then
export BAM_FILE="$alignedfiltered_bam"
source /opt/cellranger-6.1.2/sourceme.bash
samtools view -H $BAM_FILE > HFA_v5/tmp_20230711/SAM_header
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILE > HFA_v5/tmp_20230711/filtered_SAM_body
cat HFA_v5/tmp_20230711/SAM_header HFA_v5/tmp_20230711/filtered_SAM_body > HFA_v5/tmp_20230711/filtered.sam
samtools view -b HFA_v5/tmp_20230711/filtered.sam > $parsed_bam

samtools index $parsed_bam
fi


done
#STAR --genomeDir /media/NAS1/Sammi_NAS1/ref/rna-GRCh38/refdata-gex-GRCh38-2020-A/star/ --genomeLoad Remove

for FILE in HFA_v5/collapsed_bam_HFA_20230711/*.bam

do
in_file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $in_file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
out_file_name="$(echo $in_file_name | cut -d'_' -f1)_$(echo $in_file_name | cut -d'_' -f2)_$(echo $in_file_name | cut -d'_' -f3)_$cell_type"
out_file_name="$(echo $out_file_name)_tag.bam"

echo $in_file_name
echo $out_file_name

if [ ! -f "HFA_v5/tagged_bam_HFA_20230711/$out_file_name" ]; then
python HFA_tag_20230601.py -i HFA_v5/collapsed_bam_HFA_20230711/$in_file_name -o HFA_v5/tagged_bam_HFA_20230711/$out_file_name
fi

done

for FILE in HFA_v5/tagged_bam_HFA_20230711/*.bam
do

tagged_file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $tagged_file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
dedupped_file_name="$(echo $tagged_file_name | cut -d'_' -f1)_$(echo $tagged_file_name | cut -d'_' -f2)_$(echo $tagged_file_name | cut -d'_' -f3)_$cell_type"
dedupped_file_name="$(echo $dedupped_file_name)_dedup.bam"

echo $tagged_file_name
echo $dedupped_file_name


#umi_tools dedup -I collapsed_bam_HFA_20230711/H27472_Eye_SRR12355470_Astrocytes.bam -S dedupped_bam_HFA_20230711/H27472_Eye_SRR12355470_Astrocytes_dedup.bam --umi-separator=",([^;]*)" --method=unique --extract-umi-method=read_id

if [ ! -f "HFA_v5/dedupped_bam_HFA_20230711/$dedupped_file_name" ]; then
samtools index HFA_v5/tagged_bam_HFA_20230711/$tagged_file_name
/media/Home_Raid1_Voyager/sammi/.local/bin/umi_tools dedup -I HFA_v5/tagged_bam_HFA_20230711/$tagged_file_name -S HFA_v5/dedupped_bam_HFA_20230711/$dedupped_file_name --umi-tag=UB --cell-tag=CB --method=unique --extract-umi-method=tag
fi
done

for FILE in HFA_v5/dedupped_bam_HFA_20230711/*.bam

do
file_name="$(echo $FILE | cut -d'/' -f3)"

cell_type="$(echo $file_name | cut -d'_' -f4)"
cell_type="$(echo $cell_type | cut -d'.' -f1)"
new_file_name="$(echo $file_name | cut -d'_' -f1)_$(echo $file_name | cut -d'_' -f2)_$(echo $file_name | cut -d'_' -f3)_$cell_type"

aei_input_file="$(echo $new_file_name).headerfixed.bam"

echo $file_name
echo $new_file_name
echo $aei_input_file

# only need to fix header for sc
samtools view -H HFA_v5/dedupped_bam_HFA_20230711/$file_name | sed '2,195 s/SN:GRCh38_/SN:/' | samtools reheader - HFA_v5/dedupped_bam_HFA_20230711/$file_name > AEI_HFA_v7/tmp_20230711/$new_file_name.headerfixed.bam

done

/media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/RNAEditingIndex -d /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v7/tmp_20230711 \
-l /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v7/log_20230711/AEI_HFA_v7_log \
-os /media/Scratch_SSD_Voyager/sammi/RNA_editing/AEI_HFA_v7/output_index_20230711/AEI_HFA_v7_outputindex \
-f .headerfixed.bam  \
-o /media/Home_Raid1_Voyager/sammi/sammi_Scratch/RNA_editing/AEI_HFA_v7/output_mpileup_20230711/AEI_HFA_v7_outputmpileup \
--genes_expression /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz \
--refseq /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz \
--snps /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz \
-gf /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Genomes/HomoSapiens/ucscHg38Genome.fa \
-rb /media/Home_Raid1_Voyager/sammi/anaconda3/ENTER/envs/aei_v1/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/ucscHg38Alu.bed.gz \
--genome hg38 \
-mm AllMismatches \
--stranded \
--verbose
