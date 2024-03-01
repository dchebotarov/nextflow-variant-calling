#!/bin/bash

script="bam2gvcf/bam2gvcf.1rg.nf"

samplesheet="bam2gvcf/bam-sheet-test.csv"

config="nextflow.config.local"

echo Running:  $script
echo
nextflow run   \
	 -c $config  \
	 $script  \
	 -profile local_docker  \
	 -work-dir   workdir_test  \
	 --samplesheet $samplesheet 

# Reference genome can be specified in the nextflow.config file 
#  or as a command line option:
	 #   --ref references/genome1_IRGSP-1.0.fasta

