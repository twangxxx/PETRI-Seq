## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated March 2020 

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

def log_to_table(sample,bc_round):
        print(bc_round)
        if bc_round != 'bc1':
            print('only run with bc1')
            exit()
        log = open(sample + '_logs/sc_pipeline_11/' + bc_round + '.log')
        table = open(sample + '_' + bc_round + '_table.txt','w')
        table.write(bc_round + ' num' + '\t' + 'count' + '\n')
        next_line = False
        for line in log:
		if('Command line parameters' in line):
                	bcs = line[line.find(sample + '_' + bc_round + '/') + len(sample) + len(bc_round) + 2:line.find('{name}')]
                        bcs = line.split('{name}')[1].split('.fastq')[0]
                if('bc' in line):
                        line_bc = line[line.find('bc'):line.find('bc') + 6]
		if('Sequence: ' in line):
			bc_count = line.split(':')[4]
                        name = sample + '_' + line_bc + bcs
                        name = name.replace(' ','')
                        table.write(name + '\t' + bc_count.split(' ')[1] + '\n')

def log_plot(sample,bc_round,cell_num):
	#bc_plot = {}
	table = pd.read_csv(sample+'_' + bc_round + '_table.txt', sep='\t')
	#bc_plot[sample + '_' + bc_round] = table['count'].plot(kind='hist', title = sample + ' ' + bc_round + ', total = ' + str(sum(table['count'])))
	#bc_plot[sample + '_' + bc_round].get_figure().savefig(sample + '_' + bc_round + '_demultiplexing.pdf')
	#plt.cla()
	table['logcount'] = np.log10(table['count']+1)
	fig = table['logcount'].hist(bins=100,log=True,grid=False)
        #fig = table['logcount'].plot(kind='hist', title = sample + ' ' + bc_round + ', total reads = ' + str(sum(table['count'])))
        if cell_num > 0:
            if(os.path.exists(sample + '_' + bc_round + '_cumulative_frequency_table.txt')):
                table = pd.read_csv(sample + '_' + bc_round + '_cumulative_frequency_table.txt',sep='\t')
                threshold_table = table.loc[table['index'] == int(cell_num)]
                fig.axvline(np.log10(int(threshold_table['count'])), color = 'k', linestyle = '--')
        fig.set_title(sample + ' ' + bc_round + ', total reads = ' + str(sum(table['count'])))
        fig.set_xlabel('log10(number of reads)')
	fig.set_ylabel('frequency')
	fig.get_figure().savefig(sample + '_' + bc_round + '_ReadsPerBC.eps')
        plt.cla()

def freq_plot(sample,bc_round,cell_num):
	if(os.path.exists(sample + '_' + bc_round + '_cumulative_frequency_table.txt')):
		table = pd.read_csv(sample + '_' + bc_round + '_cumulative_frequency_table.txt',sep='\t') 
	else:	
		table = pd.read_csv(sample+'_' + bc_round + '_table.txt', sep='\t')
		table = table.sort_values(by='count',ascending=False)
		table['frac'] = table['count']/table['count'].sum()
		cumulative_list = np.zeros([len(table['frac'])])
		for i in range(0,len(table['frac'])):
			cumulative_list[i] = table.iloc[i]['frac'] + cumulative_list[i-1]
		table['cumulative'] = cumulative_list
		table['index'] = range(0,len(table['count']))
		table.to_csv(sample + '_' + bc_round + '_cumulative_frequency_table.txt', sep = '\t')
	threshold = table['cumulative'] < 0.9999
        fig, ax = plt.subplots(figsize=(4,4))
        ax.plot(table[threshold]['index'],table[threshold]['cumulative'],'o', rasterized=True)
        if cell_num > 0:
            ax.axvline(int(cell_num), color = 'k', linestyle = '--')
	#ax.set_title('Knee Plot, ' + sample + ' ' + bc_round)
	ax.set_xlabel('Barcode (ordered largest to smallest)',size=12)
	ax.set_ylabel('Cumulative fraction of reads',size=12)
	plt.tight_layout()
        fig.savefig(sample + '_' + bc_round + '_kneePlot.eps')
