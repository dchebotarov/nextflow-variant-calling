# nextflow-variant-calling

This is a GATK based variant calling workflow written in Nextflow.

There are two versions of the workflow:
1) a simpler version, working on samples with only one read group per sample
2) a full version for samples that may contain multiple read groups.

### Simple version

The workflow is contained in the `main.nf` script.

The workflow has the following inputs:

 - reference genome, indexed with bwa index: the path to the FASTA file (local or S3) can be given in the "--ref" commandline option
 - sample sheet CSV file with header "sample,fastq\_1,fastq\_2"

### Full version

Please see the files in the folder `multi-readgroup` 
The difference with the simple workflow is that the sample sheet need to include two additional columns:

  - column 4: read group name
  - column 5: number of read groups per sample

The last column is needed for Nextflow to group read-group BAM files as soon as the needed number is processed, in order to avoid stalling the pipeline.


 Example sample sheets are included.

 To run the workflow, please modify `run-ec2-test.sh` or `run-local-test.sh` with paths to your samplesheet and config files.
 


