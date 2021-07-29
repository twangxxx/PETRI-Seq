Created by Sydney Blattman as a supplement to "Prokaryotic single-cell RNA sequencing by in situ combinatorial indexing" (doi: 10.1038/s41564-020-0729-6)
Tavazoie Lab, Columbia University
Contact sb3292@cumc.columbia.edu with questions about code
Last updated: April 2021

Dependencies:
python 2.7.15
fastqc (v0.11.9)
cutadapt (v1.18)
GNU parallel (20190522)
UMI-tools (v0.5.5)
featureCounts (v1.6.3)
bwa (0.7.17-r1198-dirty)

Pipeline:
Example commands can be run from demo folder
1. Create a folder with the sample name (eg first1000, which is already in demo folder) and put fastq files in this folder (all lanes, paired end)
Run the following from directory containing the folder with sample name 
2. Make sure all fastq files have names in the form: {sample name}_{S#}_L00{#}_R1_001.fastq.gz or {sample name}_{S#}_L00{#}_R2_001.fastq.gz
	For example, the file names for demo sample lane 1 would be: first1000_S5_L001_R1_001.fastq.gz and first1000_S5_L001_R2_001.fastq.gz
3. SCRIPT: python [path]/PETRI_seq_scripts/scripts/sc_pipeline_11.py {full_sample} {n_lanes} 
	# full_sample is sample name and S number (eg first1000_S5)
	# n_lanes is the number of sequencing lanes for analysis - if lanes are merged, then the single file should be names with suffix _L001_R1_001.fastq.gz and n_lanes set to 1. The script will count lanes from 1 to n_lanes so always start numbering from 1.
	# sc_pipeline_11 runs fastqc, quality filter, and barcode demultiplexing
	# This step will take a few hours for ~50 million reads
3. Look at {sample}_bc1_ReadsPerBC.eps and {sample}_bc1_kneePlot.eps to determine number of BCs to include in further analysis
4. SCRIPT: [path]/PETRI_seq_scripts/scripts/pipeline.sh {sample} {n_BCs} {fasta} {gff} {custom_name} 
	# n_BCs is number of BCs to include in further analysis (typically 10000-40000)
	# fasta is location and name of fasta for alignment
	# gff is location and name of gff for feature calling - see example.gff for example format (specifically, gene names should be indicated by 'name=')
	# Custom name is a new name for the sample, corresponding to maybe the gff used or other specific input of the pipeline. For example, we might analyze the same cells by operon or by gene and would indicate that in the custom name. Custom name can be the same as sample name if desired. 
	# pipeline.sh includes a number of cleanup commands at the end. If interested in intermediate files, these can be easily commented out.

DEMO
We have provided a demo folder as an example of how to run the pipeline. We have included the first 1000 lines of 2 lanes of one of our PETRI-seq libraries. Execute the commands below to process 50 BCs from the demo library:

cd demo
python ../scripts/sc_pipeline_11.py first1000_S5 2
../scripts/pipeline.sh first1000 50 ../scripts/U00096_JE2.fa ../scripts/example.gff first1000_operons
