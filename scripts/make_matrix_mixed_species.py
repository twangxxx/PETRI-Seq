## Supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
## Written by Sydney Blattman
## Tavazoie Lab, Columbia University
## Last updated March 2020 

import sys
import pandas as pd

sample = sys.argv[1]
table = pd.read_csv(sample + '_filtered_mapped_UMIs.txt',sep='\t',index_col=0)
matrix = table[['contig:gene','UMI']]
gene_matrix = matrix.groupby(['Cell Barcode','contig:gene']).count()
gene_matrix = gene_matrix.unstack(level='Cell Barcode')
gene_matrix = gene_matrix.fillna(0)
gene_matrix = gene_matrix.transpose()
gene_matrix = gene_matrix.droplevel(0)
gene_matrix = gene_matrix.loc[:,~gene_matrix.columns.str.contains('ambiguous')]
gene_matrix.to_csv(sample + '_mixed_species_gene_matrix.txt',sep='\t')

