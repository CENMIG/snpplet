# `snpplet`

A pipeline for short variant calling from paired-end short-read genome sequencing data


## Description

Whole genome sequencing data are used for various purposes including molecular typing, classification and phylogenetic analysis. The `snpplet` pipeline processes short-read sequencing data to obtain short variants ([SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) and [indels](https://en.wikipedia.org/wiki/Indel)). It implements tools for raw read data quality control and preprocessing, read mapping, variant calling and filtering to produce an analysis-ready vcf file. The `snpplet` pipeline incorporates workflows from [samtools](http://www.htslib.org/workflow/) and [the GATK Best Practices for germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). It also includes a step to parse the SNV table into a multiple sequence alignment (MSA) in the fasta format.


## Inputs

1. Paired-end short reads (fastq format), e.g. each sample has a pair of files `<sample>_1.fastq.gz` and `<sample>_2.fastq.gz`
2. Reference genome sequence (fasta format)


## Outputs

The pipeline will generate outputs in five subdirectories under the directory `results`. Here is a brief description of output files.

Directory | Files  | Description
--------- | ------ | -----------
`multiqc_fastqc` | `before.html` | Report of the quality of raw reads of all input samples
`multiqc_fastqc` | `after.html` | Report of the quality of processed reads of all input samples after trimming using Trimmomatic
`bam` | `<sample>.bam`, `<sample>.bam.bai` | Per-sample mapped read files in BAM format and corresponding index files (`.bai`)
`bam` | `coverage.tsv`  | Summary table of mean depth and genomic coverage from `samtools coverage`
`vcf` | `<sample>.g.vcf.gz`, `<sample>.g.vcf.gz.tbi` | Per-sample vcf (in gvcf format) and index files
`vcf_joint` | `joint_filtered.vcf.gz`, `joint_filtered.vcf.gz.tbi` | Analysis-ready multi-sample vcf of high-quality SNPs
`vcf_joint` | `joint_filtered_stats.txt` | Summary statistics of the vcf file `joint_filtered.vcf.gz`
`aln` | `aln.fasta` | Multiple sequence alignment (MSA) in fasta format
`aln` | `pos.txt` | List SNP positions in the reference genome corrdinates



## Installation and usage

`snpplet` is designed and developed under Linux systems and MacOSX. For Windows users, please use a virtual machine or Windows subsystem for Linux.

1. Install dependencies listed in __Table 1__.
2. Download `snpplet`.

     ```sh
     $ git clone https://github.com/CENMIG/snpplet.git
     $ cd snpplet
     ```
3. Specify paths and parameters (see __Table 2__ below) in the `snpplet.sh` file using text editor.
    <br>Example:

    ```sh
     nextflow run snpplet.sh -resume --ref_genome 'ref/ref.fasta' --reads 'reads/*_{1,2}.fastq.gz' --results 'results' --ref_genome_name 'NC_000962.3' --stop 'true'
     ```
4. Lauch the pipeline execution.
     ```sh
     $ ./snpplet.sh
	or
	$ sh snpplet.sh
     ```
     
     <br> Paramters and running processes and its progress will be shown as figure below. For more details on how `snpplet` works, see the next section.
     ```sh
     N E X T F L O W  ~  version 20.07.1
     Launching `main.nf` [ridiculous_lorenz] - revision: f29238c695
     snpplet v0.2
     ====================================================
     Location of reference genome file: /home/bb_20/ref/ref.fasta
     Reference regenome name: NC_000962.3
     Location of raw read files (fastq): /home/bb_20/reads/*_{1,2}.fastq.gz
     Location of results: /home/bb_20/output
     Trimming option (for trimmomatic): SLIDINGWINDOW:4:30 MINLEN:70
     Mapping option (for BWA MEM): -c 100 -M -T 50
     Option for HaplotypeCaller (per-sample variant calling): -ploidy 1 -mbq 20
     Option for variant filtering (using VariantFilatration): QD < 2.0 || MQ < 40.0
     (Optional) Stop after finishing step 3 (per-sample variant calling by HaplotypeCaller): false
     =====================================================
     executor >  local (15)
     [0a/10fa43] process > PREPARE_GENOME_SAMTOOLS (ref)   [100%] 1 of 1 ✔
     [63/a74c11] process > PREPARE_GENOME_PICARD (ref)     [100%] 1 of 1 ✔
     [0c/2dbb31] process > PREPARE_GENOME_BWA (ref)        [100%] 1 of 1 ✔
     [71/d4bb14] process > FASTQC_BEFORE_TRIM (ERR1034689) [ 80%] 4 of 5
     [-        ] process > MULTIQC_FASTQC_BEFORE_TRIM      -
     [f7/9ffdf3] process > TRIM (ERR1034687)               [ 80%] 4 of 5
     [99/ee416e] process > FASTQC_AFTER_TRIM (ERR1034686)  [  0%] 0 of 4
     [-        ] process > MULTIQC_FASTQC_AFTER_TRIM       -
     [5c/7882e7] process > READ_MAPPING_BWA (ERR1034686)   [  0%] 0 of 4
     [-        ] process > COVERAGE_OUTPUT                 -
     [-        ] process > CALL_VARIANTS                   -
     [-        ] process > CREATE_SAMPLE_MAP               -
     [-        ] process > JOINT_GENOTYPING                -
     [-        ] process > FILTER_VARIANTS                 -
     [-        ] process > VCF_TO_FASTA                    -			`				

     ```
     <br>Once the run has finised, `snpplet` will report duration time that was used (example below). It will generate cache files in the directory `work/`for resuming the pipeline. You will want to delete this directory once the run has successfully completed.

    ```sh
     N E X T F L O W  ~  version 20.07.1
     Launching `main.nf` [ridiculous_lorenz] - revision: f29238c695
     snpplet v0.2
     ====================================================
     Location of reference genome file: /home/bb_20/ref/ref.fasta
     Reference regenome name: NC_000962.3
     Location of raw read files (fastq): /home/bb_20/reads/*_{1,2}.fastq.gz
     Location of results: /home/bb_20/output
     Trimming option (for trimmomatic): SLIDINGWINDOW:4:30 MINLEN:70
     Mapping option (for BWA MEM); -c 100 -M -T 50
     Option for HaplotypeCaller (per-sample variant calling); -ploidy 1 -mbq 20
     Option for variant filtering (using VariantFilatration); QD < 2.0 || MQ < 40.0
     (Optional) Stop after finishing step 4 (per-sample variant calling by HaplotypeCaller): false
     =====================================================
     executor >  local (35)
     [0a/10fa43] process > PREPARE_GENOME_SAMTOOLS (ref)   [100%] 1 of 1 ✔
     [63/a74c11] process > PREPARE_GENOME_PICARD (ref)     [100%] 1 of 1 ✔
     [0c/2dbb31] process > PREPARE_GENOME_BWA (ref)        [100%] 1 of 1 ✔
     [c5/ecda54] process > FASTQC_BEFORE_TRIM (ERR1034688) [100%] 5 of 5 ✔
     [66/1b5299] process > MULTIQC_FASTQC_BEFORE_TRIM      [100%] 1 of 1 ✔
     [9d/5fbaae] process > TRIM (ERR1034688)               [100%] 5 of 5 ✔
     [66/4a654e] process > FASTQC_AFTER_TRIM (ERR1034688)  [100%] 5 of 5 ✔
     [ba/b8909b] process > MULTIQC_FASTQC_AFTER_TRIM       [100%] 1 of 1 ✔
     [c6/d877d1] process > READ_MAPPING_BWA (ERR1034688)   [100%] 5 of 5 ✔
     [01/b30c49] process > COVERAGE_OUTPUT                 [100%] 1 of 1 ✔
     [be/aa02d7] process > CALL_VARIANTS (ERR1034685)      [100%] 5 of 5 ✔
     [1c/ebed9d] process > CREATE_SAMPLE_MAP               [100%] 1 of 1 ✔
     [b2/6117d7] process > JOINT_GENOTYPING (1)            [100%] 1 of 1 ✔
     [3b/806024] process > FILTER_VARIANTS (1)             [100%] 1 of 1 ✔
     [b0/03d1d3] process > VCF_TO_FASTA (1)                [100%] 1 of 1 ✔
     Completed at: 02-Nov-2020 12:36:48
     Duration    : 11m 11s
     CPU hours   : 0.5
     Succeeded   : 35


     Time used: 00:00:00
     ```


## Steps in the `snpplet` pipeline

#### 1. Preprocessing of NGS reads
`Trimmomatic` is used to trim Illumina adapters and low quality base calls in reads. Reads that are shorter than the specified minimum length after trimming will be dropped. The default option is `-PE -phred33, ILLUMINACLIP:TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:4:30, MINLEN:70`. This means bases from the 3' end of the read is removed if the average base quality across 4 bases is below 30, and the remaining read has at least 70 bases.

`FastQC` is used to assess the quality of the reads before and after processing with `Trimmomatic`. Then `MultiQC` is used to parse results from `FastQC` for each sample to produce summary tables of all samples.

#### 2. Mapping, indexing, sorting, and calculating depth and coverage of mapped reads
Processed reads are mapped to the reference genome sequence with `BWA MEM`. Default setting is `-c 100 -M -T 50`. The resulting SAM files are converted into BAM files using `samtools` commands (per [samtools workflow](http://www.htslib.org/workflow/)). Duplicate reads are marked by `Picard`'s `MarkDuplicates` in `GATK`. The average depth and genomic coverage of mapped reads for each sample are calculated from the BAM file using `samtools coverage`. The results are parsed into a table called `coverage.tsv`.

#### 3. Per-sample variant calling
BAM files are passed to `GATK HaplotypeCaller` for calling short variants (SNPs and indels). Default option uses the minimum base quality of 20 and ploidy of 1 for haploid organisms (`-ploidy 1 -mbq 20`). This step generates an intermediate file in the GVCF format for each sample. 

A tab-delimited file, called `sample_map.txt` containing a list of sample names and GVCF file names, is generated and saved to the _input directory_. It file looks like

            Sample_1 Sample_1.g.vcf.gz
            Sample_2 Sample_2.g.vcf.gz
            Sample_3 Sample_3.g.vcf.gz
                    .         .
            Sample_n  Sample_n.g.vcf.gz

The purpose of this file is to allow the user to manually exclude samples, for instance, based on inspection of the read quality (HTML reports from MultiQC) and mapping depth and coverage (`coverage.tsv`). For this, the user can use the option `--stop true`. Useful criteria include minimum read quality score, %GC content, minimum mean depth and minimum genome coverage. By default (`--stop false`), `snpplet` will proceed to the next step and perform joint genotyping using all samples. The list of the samples to use (i.e. not excluded) must be in the file called `sample_map_usr.txt` _in the input directory_ (e.g. created by copying from `sample_map.txt` and deleted the lines corresponding to the samples to exclude, one sample per line).

To resume the pipeline, simply execute `./sh snpplet.sh`.


#### 4. Joint genotyping
`GATK GenomicsDBImport` is used to create a GenomicsDB workspace and then import single-sample GVCF files into this workspace. Then `GATK GenotypeGVCFs` is used to perform joint genotyping on the GenomicsDB workspace to produce a single multi-sample VCF file. By default, `GenomicsDBImport` uses all samples as listed in the `sample_map.txt`. If the user-created file `sample_map_usr.txt` is present in the input directory, the list of samples in this file will be used instead.

										
#### 5. Variant filtering
By default, indels are filtered out using `GATK SelectVariants` with option `-select-type SNP`. The low-quality SNPs are also filtered out using `GATK VariantsFiltration` with option `--exclude-filtered`. The default filter is
`QD < 2.0 || MQ < 40.0`. After filtering, the resulting analysis-ready VCF file is generated as `joint_filtered.vcf.gz` in the `vcf_joint` directory in output directory.


#### 6. VCF to FASTA
Finally, the filtered VCF file `joint_filtered.vcf.gz` is processed to produce a multiple sequence alignment file called `aln.fasta` (FASTA format, in `aln` directory in the output directory) using a combination of tools including `GATK VariantsToTable` and `datamash`. Positions of the SNPs in the reference coordiates are reported a separate file named  `pos.txt` in the same directory.


### Table 1 - Dependencies

* [Nextflow](https://www.nextflow.io/) (v20.07 or later)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [BWA](http://bio-bwa.sourceforge.net/)
* [Samtools](http://www.htslib.org/)
* [BCFtools](http://www.htslib.org/)
* [GATK](https://gatk.broadinstitute.org/hc/en-us) (v4.1.8.1 or later)
* [Datamash](https://www.gnu.org/software/datamash/)
* Java for Nextflow and GATK. Please check specific requirements for Nextflow and GATK.


### Table 2 - Arguments (paths and parameters)

Argument | Description   | Default value
-------- | ------------- | -------------
`--reads` | the location of the read fastq files | `reads/*_{1,2}.fastq.gz`
`--ref_genome` | the location of the read fastq files | `ref/reference.fasta`
`--results` | the location where the results will be stored |  `results`
`--ref_genome_name` | genome name | `NC000962.3`
`--trimming_option` | parameter for trimming (Trimmomatic) | `"SLIDINGWINDOW:4:30 MINLEN:70"`
`--mapping_option` | parameter for read mapping (BWA) | `"-c 100 -M -T 50"`
`--haplotypecaller_option` | parameter for per-sample variant calling (GATK HaplotypeCaller) | `"-ploidy 1 -mbq 20"`
`--genomicsdbimport_option` | parameter for joint (GATK GenomicsDBImport) | `"--batch-size 200"`
`--genotypegvcfs_option` | parameter for joint (GATK GenotypeGVCFs) | `"-ploidy 1"`
`--variant_filter` | parameter for variant filtering (GATK VariantsFiltration) |  `"QD < 2.0 || MQ < 40.0"`
`--variant_filter_name` | name of the variant filter | `"qd-mq"`
`--selectvariant_option` | parameter for selecting variants (GATK SelectVariants) | `"--exclude-filtered -select-type SNP"`
`--stop` | whether to pause after step 3 (per-sample variant calling) | `false`


## Quick installation of `snpplet`'s dependencies

For Linux systems using bash shell.

#### 1. Install `homebrew` and `homebrew/science/bio`

1.1 Install `curl` and `git`:

          
          $ sudo apt-get update
          

1.2 Install [`homebrew`](https://brew.sh/):

         
         $ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)" 
         

1.2.1 Set `$PATH`:

         
         $ echo 'eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv)' >> /home/user/.profile 
         

1.2.2 Evaluate `.profile`:

         
         $ eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv) 
         

1.3 Install `build-essential`:

          
          $ sudo apt-get install build-essential 
          

1.4 Install `gcc`:

        
          $ brew install gcc 
          

1.5 Tap `brewsci/bio`:

          
          $ brew tap brewsci/bio 
          

#### 2. Install dependencies for `snpplet` via `homebrew`

     
     $ brew install nextflow gatk fastqc samtools bcftools trimmomatic bwa datamash
     

Launch `nexflow` to install its dependencies by running `nextflow`


#### 3. Install `multiqc` (a python package) using `pip`:

        
        $ pip3 install multiqc
        



# Schematic outline
![workflow](https://github.com/CENMIG/snpplet/blob/master/figures/steps.png)

You can cite the snpplet zenodo record for a specific release using the following <a href="https://zenodo.org/badge/latestdoi/310209417"><img src="https://zenodo.org/badge/310209417.svg" alt="DOI"></a>

