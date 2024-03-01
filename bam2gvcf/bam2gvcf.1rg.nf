// Variant calling pipeline part
// 
//  from sorted BAM to gVCF
// (no gVCF merging)
// Assume each sample has exactly one BAM file (no multiple read groups per sample)
// The sample names and BAM files must be given in  samplesheet.csv file 
//  with header line "sample,bam"

// Using Nextflow DSL 2

//   Data parameters: sample sheet and reference files basename
// can be set in the nextflow.config file 
// or in the commandline arguments

// params.samplesheet = 'samplesheet.csv'

ref = params.ref

println "ref:  $ref"  

params.ref_fai = "${params.ref}.fai"

//refBase = ref.baseName
refBase = ref.take(ref.lastIndexOf('.'))
params.ref_dict = file(refBase + '.dict')

// Additional parameters for read group setting
params.CN = 'TXG' // sequencing center
params.PL = "ILLUMINA" // sequencing platform

// ########### SOFTWARE AND CONTAINERS ############
// gatk_docker = "public.ecr.aws/biocontainers/gatk4:4.1.9.0--py39_0"
gatk_invoc = params.gatk_invoc          // "gatk"

// picard_docker = gatk_docker // since gatk v4
picard_invoc = params.picard_invoc     // "java -jar /usr/gitc/picard.jar "

bwa_invoc = params.bwa_invoc           //    "/usr/gitc/bwa "
samtools_invoc = params.samtools_invoc // "samtools"


// ############################################################### //
// ##############   Prepare inputs ############################### //
// ############################################################### //
// Parsing sample sheet, sending to channel "samples_ch"
Channel.fromPath(params.samplesheet)
  .splitCsv(header: true)
  .map { row -> tuple( row.sample, file(row.bam) ) }
  .set { samples_ch }

// channels for reference 
ref_rel_files = Channel.value( params.ref )
            .map { f -> tuple( 
                 file(f), 
                 file("${f}.amb"), 
                 file("${f}.ann"), 
                 file("${f}.bwt"), 
                 file("${f}.pac"),
                 file("${f}.sa")
                 )  }

ref_files_for_varcall = Channel.value( tuple(
               file(params.ref), 
               file(params.ref_fai), 
               file(params.ref_dict)
               ) )

// ####################################################### //
//                 PROCESS - ALIGN (not used)
// ####################################################### //
process align {
 input:
 tuple path(ref), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa)
 tuple val(sid), path(read1), path(read2)
 output:
 tuple val(sid), path("${sid}.bam"), emit: bam
 script:
 """
 echo "cpu: ${task.cpus}  mem: ${task.memory}"
 echo "${workflow.containerEngine}"
 echo "##   "  ${sid} ':' ${read1} ${read2}
 $bwa_invoc  mem -t ${task.cpus}  ${ref} ${read1} ${read2} > out.sam
 samtools view -Obam out.sam > ${sid}.bam
 """
}


// ####################################################### //
// ##### Processing BAMs 
// ####################################################### //

process bamProcess {
 time '12h'

 publishDir 'stats', pattern: "*.metrics", mode: 'move'
 publishDir 'final-bams', mode: 'copy', pattern: "*.proc.ba*"
 // storeDir 'processed-bams'  // mode: 'copy'

 // container "${picard_docker}" // container settings in config file

 input:
 tuple val(sid), path(bam)
 
 output:
 tuple val(sid), path("${sid}*.proc.bam"), path("${sid}*.proc.bai"),  emit: proc_bam
 path "*.metrics",  emit: mkdup_metrics_file, optional: true
 //path("${sid}*.bai"), optional: true

 script:
 """
 echo $bam 
  # sort sam
  $picard_invoc  SortSam -I $bam -O sorted.bam \
    --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    --SORT_ORDER coordinate 
  # fix mate
  $picard_invoc  FixMateInformation -I  sorted.bam  -O fxmt.bam \
      --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE

  # mark dup
  $picard_invoc  MarkDuplicates -I  fxmt.bam  -O markdup.bam \
      --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    --METRICS_FILE ${sid}.markdup.metrics \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000

  # add/replace read groups
  $picard_invoc  AddOrReplaceReadGroups -I  markdup.bam  -O addrep.bam \
      --VALIDATION_STRINGENCY LENIENT \
    --CREATE_INDEX TRUE \
    --RGLB lib1 \
    --RGSM ${sid} \
    --RGID ${sid} \
    --RGPU PU1 \
    --RGPL ${params.PL} \
    --RGCN ${params.CN}

  # set the name of final file
  mv addrep.bam ${sid}.proc.bam
  mv addrep.bai ${sid}.proc.bai
 """

}



// ######################################################### //
// ######        Variant Calling  ##
// ######################################################## //

process variantCall {
 time '2d'

 publishDir "${params.outdir}/gvcf/", pattern: "*vcf*", mode: 'copy'
 // storeDir 'store-gvcf'
 publishDir "${params.outdir}/bamout/", pattern: "*bamOut.ba*", mode: 'copy'

 //  publishDir 'results-gvcf', mode: 'copy'
 // storeDir 'store-gvcf'

 // container "$gatk_docker"
 // containerOptions ' -v `pwd`:/data -v results-gvcf:/gatk/results-gvcf '

 input:
 tuple path(ref), path(ref_fai), path(ref_dict)
 tuple val(sid), path(bam), path(bai)
 

 output:
 tuple val(sid), path("*.vcf.gz"), emit: gvcf
 path("*.bamOut.ba*"),  emit: bamOut

 script:

  def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

 """
 ## for older GATK version
 ##optionalT=" -T "
 
 bands="--gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30  --gvcf-gq-bands 33 --gvcf-gq-bands 35 --gvcf-gq-bands 37  --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 42  --gvcf-gq-bands 45 --gvcf-gq-bands 46  --gvcf-gq-bands 48  --gvcf-gq-bands 50 --gvcf-gq-bands 52  --gvcf-gq-bands 55  --gvcf-gq-bands 57  --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99"

 # consider adding --gvcf-bands
 

 echo "#----------- Working directory ------------#"
 pwd
 echo "#------------------------------------------#"

 $gatk_invoc --java-options "-Xmx${avail_mem}M -XX:-UsePerfData"  HaplotypeCaller \
         -R $ref \
         -ERC  GVCF \
         -I $bam \
         -O ${sid}.g.vcf.gz  \
         -mbq 20 \
          --showHidden true \
         -bamout ${sid}.bamOut.bam \
  
  # 
  #      \$bands

 """

}


// MAIN WORKFLOW 
workflow {
  
  main:
    
    // align( ref_rel_files, samples_ch)   
    
    // bamProcess(align.out.bam)

    bamProcess(samples_ch)
    variantCall( ref_files_for_varcall, bamProcess.out.proc_bam )

    // align.out.verbiage.view()

 emit:
   variantCall.out.gvcf
}
