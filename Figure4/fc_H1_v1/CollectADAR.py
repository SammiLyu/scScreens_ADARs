### adapted from Adriano Fonzino et al.
## (ENSG00000160710) ADAR
## (ENSG00000197381) ADARB1
## (ENSG00000185736) ADARB2

import sys, os
print('Sample\tGene\tTPM')
f=open('/media/Scratch_SSD_Voyager/sammi/RNA_editing/fc_v1/counts_all_tpm.txt')
tpm_list={}
for i in f: 
  if i.startswith('Geneid'):
    h=i.strip().split('\t')
    for j in range(len(h)):
      if h[j].startswith('X.media.Scratch_SSD_Voyager.sammi.RNA_editing.dedupped_bam_H1_20230102'):
        sn_in = h[j].strip().split("_")
        for k in range(len(sn_in)):
          if sn_in[k].startswith('20230102'): sn_p1 = sn_in[k]
          if sn_in[k].startswith('ter'): sn_p2 = sn_in[k]
        sn = sn_p1.replace('.', '_')+'_'+sn_p2
        tpm_list[sn] = []
        if h[j].startswith('X.media.Scratch_SSD_Voyager.sammi.RNA_editing.dedupped_bam_H1_'+sn_p1+'_'+sn_p2): tpm_list[sn].append(j)
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
