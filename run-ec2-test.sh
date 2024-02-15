

script="main.nf"

# samplesheet="sampleSheet.ec2-test.csv"  # real samples
samplesheet="sampleSheet.aws.csv"  # small region


# ref="references/genome1_IRGSP-1.0.fasta" # if local copy
ref="s3://irri-nextflow-awsbatch-test/references/genome1_IRGSP-1.0.fasta"

echo Running:  $script
echo
nextflow run   -c nextflow.config.main  \
	 $script  \
	 -profile aws_poc  \
	 -with-docker  \
	 -work-dir   workdir  \
	 --samplesheet $samplesheet \
    --ref $ref

