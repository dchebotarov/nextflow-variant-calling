

script="main.nf"

samplesheet="sampleSheet.test.csv"


echo Running:  $script
echo
nextflow run   -c nextflow.config.test  \
	 $script  \
	 -profile irrimac  \
	 -with-docker  \
	 -work-dir   workdir  \
	 --samplesheet $samplesheet \
   --ref test-data/ref/Chr7.fasta   

# --ref references/genome1_IRGSP-1.0.fasta

