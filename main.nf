testing
/*
 * Copyright (c) 2020
 * Center for Microbial Genomics, Mahidol University
 * 
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */


/* 
 * 'snpplet' - A Nextflow pipeline for variant calling from NGS data
 * 
 * Yuttapong Thawornwattana
 * Bharkbhoom Jaemsai
 */


/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */
// paths & inputs
baseDir          = "$HOME"
params.ref_genome    = "$baseDir/ref/reference.fasta"
params.reads     = "$baseDir/reads/*_{1,2}.fastq.gz"
params.results   = "$baseDir/results"
params.adapter   = "$baseDir/data/adapters"
params.ref_genome_name = "NC_000962.3"

// parameters for trimming
params.trimming_option = "SLIDINGWINDOW:4:30 MINLEN:70"

// parameters for read mapping
params.mapping_option = "-c 100 -M -T 50"

// parameters for variant calling
params.haplotypecaller_option = "-ploidy 1 -mbq 20"
params.genomicsdbimport_option = "--batch-size 200"
params.genotypegvcfs_option = "-ploidy 1"
params.sample_map_usr = "$baseDir/data/sample_map_usr.txt"  // OPTIONAL

// parameters for variant filtering
params.variant_filter="QD < 2.0 || MQ < 40.0"
params.variant_filter_name="qd-mq"
params.selectvariant_option = "--exclude-filtered -select-type SNP"

// Optional stop for sample QC before joint genotyping (Step 4 below)
params.stop = false


log.info """\
snpplet v0.2
========================================================
Reference genome sequence: $params.ref_genome
Reference genome name: $params.ref_genome_name
Raw reads: $params.reads
Results  : $params.results
Trimming option (for trimmomatic): $params.trimming_option
Mapping option (for BWA MEM); $params.mapping_option
Option for HaplotypeCaller (per-sample variant calling): $params.haplotypecaller_option
Option for variant filtering (using VariantFilatration): $params.variant_filter
(Optional) Stop after finishing step 3 (per-sample variant calling by HaplotypeCaller): $params.stop
========================================================
"""


/* 
 * Import modules 
 */
include { 
  PREPARE_GENOME_SAMTOOLS;
  PREPARE_GENOME_BWA;
  PREPARE_GENOME_PICARD;
  FASTQC_BEFORE_TRIM;
  MULTIQC_FASTQC_BEFORE_TRIM;
  TRIM;
  FASTQC_AFTER_TRIM;
  MULTIQC_FASTQC_AFTER_TRIM;
  READ_MAPPING_BWA;
  COVERAGE_OUTPUT;
  CALL_VARIANTS;
  CREATE_SAMPLE_MAP;
  JOINT_GENOTYPING;
  FILTER_VARIANTS;
  VCF_TO_FASTA
} from './modules.nf' 


/* 
 * Main pipeline steps
 */
workflow {
  // input: paired-end reads
  read_pairs = Channel.fromFilePairs(params.reads)

  // optional input: user-specified list of samples to be included
  // in joint genotyping [default: use all samples]
  sample_map_usr = Channel.fromPath(params.sample_map_usr)

  // STEP 0: Data preparation
  PREPARE_GENOME_SAMTOOLS(params.ref_genome)
  PREPARE_GENOME_PICARD(params.ref_genome)
  PREPARE_GENOME_BWA(params.ref_genome)

  // STEP 1: Quality trimming, with read QC before and after
  FASTQC_BEFORE_TRIM(read_pairs)
  MULTIQC_FASTQC_BEFORE_TRIM(FASTQC_BEFORE_TRIM.out.collect())
  TRIM(
    read_pairs,
    params.adapter,
    params.trimming_option
  )
  FASTQC_AFTER_TRIM(TRIM.out)
  MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.collect())

  // STEP 2: Read mapping using BWA MEM
  READ_MAPPING_BWA(
    params.ref_genome,
    params.ref_genome_name,
    PREPARE_GENOME_BWA.out,
    TRIM.out,
    params.mapping_option)

  // Report per-sample depth and coverage statistics
  COVERAGE_OUTPUT(
    READ_MAPPING_BWA.out[1].collect())

  // STEP 3: Per-sample variant calling using GATK HaplotypeCaller
  CALL_VARIANTS(
    params.ref_genome, 
    PREPARE_GENOME_SAMTOOLS.out,
    PREPARE_GENOME_PICARD.out, 
    READ_MAPPING_BWA.out[0].groupTuple(),
    params.haplotypecaller_option)

  // Create a list of samples to allow user to exclude samples
  // before joint genotyping
  CREATE_SAMPLE_MAP(
    CALL_VARIANTS.out.vcf.collect())

  // STEP 4: Joint genotyping using GATK GenotypeGVCFs
  JOINT_GENOTYPING( 
    params.ref_genome,
    PREPARE_GENOME_SAMTOOLS.out,
    PREPARE_GENOME_PICARD.out, 
    params.ref_genome_name, 
    CALL_VARIANTS.out.vcf.collect(),
    CALL_VARIANTS.out.vcf_tbi.collect(),
    sample_map_usr,
    params.genomicsdbimport_option,
    params.genotypegvcfs_option)

  // STEP 5: Variant filtering
  FILTER_VARIANTS(
    params.ref_genome,
    PREPARE_GENOME_SAMTOOLS.out,
    PREPARE_GENOME_PICARD.out, 
    params.variant_filter,
    params.variant_filter_name,
    JOINT_GENOTYPING.out,
    params.selectvariant_option)

  // STEP 6: Generating multiple sequence alignment from vcf
  VCF_TO_FASTA(FILTER_VARIANTS.out[0])
}
