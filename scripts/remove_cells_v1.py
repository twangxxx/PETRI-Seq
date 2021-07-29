## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated March 2020 

import pandas as pd
import sys
import os

sample = sys.argv[1] 
threshold = sys.argv[2] 

table = open(sample + '_bc1_cumulative_frequency_table.txt')
output = open(sample + '_selected_cumulative_frequency_table.txt', 'w')
table.next()

for line in table:
    name = line.split('\t')[1]
    R1_file_name = sample + '_bc1/' + sample + '_R1' + name[name.find('_bc1_'):len(name)] + '.fastq.gz'
    R1_file_name = R1_file_name.replace(' ' , '')
    R2_file_name = sample + '_bc1/' + sample + '_R2' + name[name.find('_bc1_'):len(name)] + '.fastq.gz'
    R2_file_name = R2_file_name.replace(' ' , '')
    if (int(line.split('\t')[5]) < int(threshold)):
	output.write(line)
    elif int(line.split('\t')[2]) > 0:
        os.system('rm ' + R1_file_name)
        os.system('rm ' + R2_file_name)
os.system('rm ' + sample + '_bc1/*_unknown*')
os.system('mv ' + sample + '_bc1/ ' + sample + '_selected_cells')
output.close()
table.close()
