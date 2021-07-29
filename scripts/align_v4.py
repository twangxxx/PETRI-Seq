## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated March 2020 

import sys
import os

ID = sys.argv[1]
if len(sys.argv) > 3:
    sample = sys.argv[3]
else:
    sample = ID

template = sys.argv[2]
table = open(ID+'_selected_cumulative_frequency_table.txt')

os.system('mkdir ' + sample + '_bwa_sai')
os.system('mkdir ' + sample + '_bwa_sam')
R2_list = {}
R2_list[0] = ''
i=0
group=0
for line in table:
    i+=1
    name = line.split('\t')[1]
    R2_file_name = 'R2' + name[name.find('_bc1_'):len(name)]
    R2_file_name = R2_file_name.replace(' ' , '')
    #print(R2_file_name)
    if R2_list[group] == '':
        R2_list[group] = R2_file_name
    else:
        R2_list[group] = R2_list[group] + '\n' + R2_file_name
    if  i%4000 == 0:
        group +=1
        R2_list[group] = ''

for group in R2_list:
    os.system('echo "' + R2_list[group] + '" | time parallel --bar --files --results ' + sample + '_bwa_sai -j10 bwa aln -n 1 ' + template + ' ' + sample + '_R2_trimmed/' + ID + '_{}_R2_trimmed.fastq.gz')
    os.system('echo "' + R2_list[group] + '" | time parallel --bar --files --results ' + sample + '_bwa_sam -j10 bwa samse -n 14 ' + template + ' ' + sample + '_bwa_sai/1/{}/stdout ' + sample + '_R2_trimmed/' + ID + '_{}_R2_trimmed.fastq.gz')

os.system('rm -r ' + sample + '_bwa_sai/')
