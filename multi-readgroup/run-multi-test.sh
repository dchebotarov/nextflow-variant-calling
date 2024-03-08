

script="multi-readgroup/multirg.nf"

# samplesheet="sampleSheet.multirg.csv"
samplesheet="multi-readgroup/sampleSheet.multirg.ct.csv"

tmpdir="./TMPDIR"
mkdir -p "$tmpdir"


echo Running:  $script
echo
nextflow run   -c multi-readgroup/nextflow.config  \
	 $script  \
	 -profile aws_poc  \
	 -work-dir   workdir  \
	 --samplesheet $samplesheet \
   -with-docker  \
	 --tmpdir $tmpdir \
   --ref test-data/ref/Chr7.fasta   "$@" 


