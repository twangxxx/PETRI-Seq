## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated March 2020 

import os

def uniform_sam(sample):

	input_sam = open(sample)
        #input_sam = open(sample + '.sam')
	output = open(sample+'_uniform_XT.sam','w')

	for line in input_sam:
		if("XT:A:U" in line):
			output.write(line.replace("XT:A:U", "XT:A:N"))
		elif("XT:A:R" in line):
			output.write(line.replace("XT:A:R", "XT:A:N"))
		else:
			output.write(line)
	output.close()


def no_xt(sample):

	input_sam = open(sample + '.sam')
	output = open(sample+'_no_xt.sam','w')
	for line in input_sam:
        	if("XT:" in line):
                	output.write(line.replace("XT:", "XN:"))
        	else:
                	output.write(line)
	output.close()

def no_xt_new(sample_in,sample_out):
        input_sam = open(sample_in)
        output = open(sample_out,'w')
        for line in input_sam:
                if("XT:" in line):
                        output.write(line.replace("XT:", "XN:"))
                else:
                        output.write(line)
        output.close()
