### adapted from Adriano Fonzino et al.
## (ENSG00000160710) ADAR
## (ENSG00000197381) ADARB1
## (ENSG00000185736) ADARB2

import sys, os
print('Sample\tGene\tTPM')
f=open('/media/Scratch_SSD_Voyager/sammi/RNA_editing/fc_TS_v1/counts_all_organize_tpm.txt')
tpm_list={}
for i in f: 
  if i.startswith('Geneid'):
    h=i.strip().split('\t')
    for j in range(len(h)):
      if h[j].startswith('X.media.Scratch_SSD_Voyager.sammi.RNA_editing.dedupped_bam_TS_20230125'):
        sn_in1 = h[j].strip().split("_")
        for k in range(len(sn_in1)):
          if sn_in1[k].startswith('organize'): sn_p1 = sn_in1[k]
        sn_in2 = h[j].strip().split("_wk")
        for k in range(len(sn_in2)):
          if sn_in2[k].endswith('.bam'): sn_p2 = sn_in2[k]
        sn = sn_p1.replace('.', '_')+'_wk'+sn_p2
        tpm_list[sn] = []
        if h[j].startswith('X.media.Scratch_SSD_Voyager.sammi.RNA_editing.dedupped_bam_TS_20230125_'+sn_p1+'_wk'+sn_p2): tpm_list[sn].append(j)
    continue
  l=i.strip().split('\t')
  if l[0]=="GRCh38_ENSG00000160710": 
    for j in range(len(l)):
      for k in list(tpm_list):
        if j in tpm_list[k]: print('%s\t%s\t%s' %(k,'ADAR',l[j]))
  if l[0]=="GRCh38_ENSG00000197381": 
    for j in range(len(l)):
      for k in list(tpm_list):
        if j in tpm_list[k]: print('%s\t%s\t%s' %(k,'ADARB1',l[j]))
  if l[0]=="GRCh38_ENSG00000185736": 
    for j in range(len(l)):
      for k in list(tpm_list):
        if j in tpm_list[k]: print('%s\t%s\t%s' %(k,'ADARB2',l[j]))
f.close()
