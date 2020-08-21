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
 * 'tuber' - A Nextflow pipeline for variant calling from NGS data
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
baseDir          = "$HOME/proj/pipeline"
params.genome    = "$baseDir/data/reference/ref.fasta"
params.reads     = "$baseDir/data/reads/*_{1,2}.fastq.gz"
params.results   = "$baseDir/results"
params.adapter   = "$baseDir/data/adapters"
params.genome_name = "NC_000962.3"

// parameters for trimming
params.trim_option = "SLIDINGWINDOW:4:30 MINLEN:70"

// parameters for read mapping
params.bwa_option = "-c 100 -M -T 50"

// parameters for variant calling
params.haplotypecaller_option = "-ploidy 1 -mbq 20"
params.genomicsdbimport_option = "--batch-size 200"
params.genotypegvcfs_option = "-ploidy 1"
params.sample_map_usr = "$baseDir/data/sample_map_usr.txt"

// parameters for variant filtering
params.variant_filter="QD < 2.0 || MQ < 40.0"
params.variant_filter_name="qd-mq"
params.mask_positions = "$baseDir/data/snp-mask-reduced.list"
params.selectvariant_option = "--exclude-filtered -select-type SNP"

// params.stop = true
params.stop = false


log.info """\
tuber v0.1
================================
genome   : $params.genome
reads    : $params.reads
results  : $params.results
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

  // STEP 1: Data preparation
  PREPARE_GENOME_SAMTOOLS(params.genome)
  PREPARE_GENOME_PICARD(params.genome)
  PREPARE_GENOME_BWA(params.genome)

  // STEP 2: Quality trimming, with read QC before and after
  FASTQC_BEFORE_TRIM(read_pairs)
  MULTIQC_FASTQC_BEFORE_TRIM(FASTQC_BEFORE_TRIM.out.collect())
  TRIM(
    read_pairs,
    params.adapter,
    params.trim_option
  )
  FASTQC_AFTER_TRIM(TRIM.out)
  MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.collect())

  // STEP 3: Read mapping using BWA MEM
  READ_MAPPING_BWA(
    params.genome,
    params.genome_name,
    PREPARE_GENOME_BWA.out,
    TRIM.out,
    params.bwa_option)

  // Report per-sample depth and coverage statistics
  COVERAGE_OUTPUT(
    READ_MAPPING_BWA.out[1].collect())

  // STEP 4: Variant calling
  // Per-sample variant calling using GATK HaplotypeCaller
  CALL_VARIANTS(
    params.genome, 
    PREPARE_GENOME_SAMTOOLS.out,
    PREPARE_GENOME_PICARD.out, 
    READ_MAPPING_BWA.out[0].groupTuple(),
    params.haplotypecaller_option)

  // Create a list of samples to allow user to exclude samples
  // before joint genotyping
  CREATE_SAMPLE_MAP(
    CALL_VARIANTS.out.vcf.collect())

  // Joint genotyping using GATK GenotypeGVCFs
  JOINT_GENOTYPING( 
    params.genome,
    PREPARE_GENOME_SAMTOOLS.out,
    PREPARE_GENOME_PICARD.out, 
    params.genome_name, 
    CALL_VARIANTS.out.vcf.collect(),
    CALL_VARIANTS.out.vcf_tbi.collect(),
    sample_map_usr,
    params.genomicsdbimport_option,
    params.genotypegvcfs_option)

  // Variant filtering
  FILTER_VARIANTS(
    params.genome,
    PREPARE_GENOME_SAMTOOLS.out,
    PREPARE_GENOME_PICARD.out, 
    params.variant_filter,
    params.variant_filter_name,
    params.mask_positions,
    JOINT_GENOTYPING.out,
    params.selectvariant_option)

  // STEP 5: Generating multiple sequence alignment from vcf
  VCF_TO_FASTA(FILTER_VARIANTS.out[0])
}
