

script="multi-readgroup/multirg.nf"

# samplesheet="sampleSheet.multirg.csv"
samplesheet="multi-readgroup/sampleSheet.multirg.ct.csv"

echo Running:  $script
echo
nextflow run   -c multi-readgroup/nextflow.config  \
	 $script  \
	 -profile aws_poc  \
	 -work-dir   workdir  \
	 --samplesheet $samplesheet \
   -with-docker  
   --ref test-data/ref/Chr7.fasta   "$@" 


