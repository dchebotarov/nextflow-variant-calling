# nextflow-variant-calling

This is a GATK based variant calling workflow written in Nextflow.

The workflow is contained in the `main.nf` script.

The workflow has the following inputs:

 - reference genome, indexed with bwa index: the path to the FASTA file (local or S3) can be given in the "--ref" commandline option
 - sample sheet CSV file with header "sample,fastq\_1,fastq\_2"

 Example sample sheets are included.

 To run the workflow, please modify `run-ec2-test.sh` or `run-local-test.sh` with paths to your data

