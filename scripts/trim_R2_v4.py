## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated March 2020 

import sys
import os
import multiprocessing
import pandas as pd
script_dir = sys.argv[0].split('trim_R2_v4.py')[0]
def trim_R2(sample,ID,name,bc_to_seq):
    R2_file_name = ID + '_R2' + name[name.find('_bc1_'):len(name)]
    R2_file_name = R2_file_name.replace(' ' , '')
    R2_short_name = R2_file_name.replace(ID+'_','')
    barcodes = R2_file_name.split('_bc')
    bc1 = barcodes[1].split('_')[1]
    sequence = bc_to_seq[bc1]
    os.system('cutadapt -a "' + sequence + 'TCTGGCGTAGGAGG;min_overlap=4" -a "TCTGGCGTAGGAGG;min_overlap=8" --minimum-length 16 -o ' + sample + '_R2_trimmed/' + R2_file_name + '_R2_trimmed.fastq.gz ' + sample + '_selected_cells/' + R2_file_name + '.fastq.gz >> ' + sample + '_logs/trim_R2_v4/trim_R2_v4.log')
    print(R2_file_name)

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]
def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)
def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq

sample = sys.argv[1] 
if len(sys.argv) > 2:
    ID = sys.argv[2]
else:
    ID = sample
table = open(sample+'_selected_cumulative_frequency_table.txt')
os.system('mkdir ' + sample + '_logs/trim_R2_v4')
bc1_list = open(script_dir + '/sc_barcodes_v2/BC1.fa')
bc_to_seq = {}
for line in bc1_list:
    if '>' in line:
        barcode = line.split('bc1_')[1].split('\n')[0]
    else:
        seq = reverse_complement(line.split('N$')[0])
        bc_to_seq[barcode] = seq
os.system('mkdir ' + sample + '_R2_trimmed')
for line in table:
    name = line.split('\t')[1]
    trim_R2(sample,ID,name,bc_to_seq)
