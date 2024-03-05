

script="multi-readgroup/multirg.nf"

# samplesheet="sampleSheet.multirg.csv"
samplesheet="multi-readgroup/sampleSheet.multirg.ct.csv"

echo Running:  $script
echo
nextflow run   -c nextflow.config.local  \
	 $script  \
	 -profile local_mac  \
	 -work-dir   workdir  \
	 --samplesheet $samplesheet \
   --ref test-data/ref/Chr7.fasta   "$@" 

# --ref references/genome1_IRGSP-1.0.fasta

