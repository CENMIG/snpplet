# `snpplet`
A nexflow pipeline for short variant calling from NGS data
## Description
Whole genome sequencing data are used for various purposes including molecular typing, classification and phylogenetic analysis. The `snpplet` pipeline processes NGS data to obtain short variants (SNVs and INDELs). This pipeline implements tools for raw data quality control and pre-processing, read mapping, variant calling and filtering to produce a ready to use vcf and multiple sequence alignment files. The snpplet pipeline includes workflow of GATK Best Practices for germline short variant discovery, [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-), which allows users to work on large datasets. In addition to the GATK Best Practices, snpplet includes a step to parse the SNV table into a multiple sequence alignment in fasta format. (`snpplet` runs well under Linux OS and MacOSX; for Windows users `snpplet` also runs well under linux in a virtual machine or Windows subsystem for linux)

## Inputs
<br>The pipeline needs the following files as the inputs
<br>1. Raw reads; file name must match the following scheme: _{1,2}.fastq.gz 
<br>2. Reference genome sequence (fasta format)
<br>3. (optional) Interval lists (BED format); formatted as 3-column format


## Outputs 
<br>The pipeline will produce outputs into 5 separate directories including aln, bam, multiqc_fastqc, vcf and vcf_joint under the directory /user_defined_path/results. Here is a brief description of output files:

File | Description
------ | ------ 
before.html | quality control check before trimming- html summary 
after.html | results quality control check after trimming using Trimmomatic with defined parameters 
bam | binary alignment map - binary version of SAM (sequence alignment/map)
coverage.tsv  | summary table of mean depth and genome coverage (numreads, covbases, coverage, meandepth, meanbaseq) 
*.g.vcf.gz | per-sample vcf 
joint_filtered.vcf | joint-called vcf filtered with default or user defined filter options
aln.fasta | multiple sequence alignment in fasta format (concatenated snp alignment) 
pos.txt | position of polymorphic sites corresponding positions in the reference genome sequence

## Installation and usage
1. Install dependencies listed in `table 1`.
2. Download the `snpplet` pipeline scripts.
     ```sh
     $ git clone https://github.com/CENMIG/snpplet.git
     $ cd snpplet
     ```
3. Specify paths and parameters for `nextflow` (listed in `table 2`) in the `snpplet.sh` file using text editor.
    <br>Example:
    ```sh
     nextflow run snpplet.sh -resume --ref '/home/bb/ref/ref.fasta' --reads '/home/bb/reads/*_{1,2}.fastq.gz' --results /home/bb/output_here/ --ref_genome 'NC_000962.3' --stop 'true'
     ```
4. Lauch the pipeline execution with the following command:
     ```sh
     $ ./snpplet.sh
	or
	$ sh snpplet.sh
     ```
     
     <br> Paramters that you have specified, and and processes and its progress will be shown as figure below. For more details on how `snpplet` works, please find next section.
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
     <br>Once the run has finised, `snpplet` will report duration time that was used (Example shown below). The pipeline will generate cache files in directory `/snpplet/work/`, you can delete this directory.

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
## The snpplet pipeline will perform the followings:

#### 1. Preprocessing raw reads
`Trimmomatic` will be used for clipping Illumina adapters and trimming low quality (sequences in the) reads. Reads that are shorter than the defined length will be dropped. By default it is set as: -PE -phred33, ILLUMINACLIP:TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:4:30, MINLEN:70
`Fastqc` will be used for QC checking sequencing reads before and after processing with ``Trimmomatic. Then MultiQC will parse summary stats from fastqc results of each sample to produce summary tables (multiqc_fasqc.txt and multiqc_general_stats.txt and HTML reports in multiqc_fastqc/ in output directory.

#### 2. Mapping, indexing, sorting, and calculating depth and coverage
Processed reads will be aligned to the reference genome sequence with `BWA MEM`. Default is set as: -c 100 -M -T 50. SAM files will be converted into BAM files by `Samtools` and duplicate reads will be marked by `Picard`. Average depth and genome coverage of each sample are calculated from BAM using Samtools, results will be parsed into coverage.tsv file in bam/ in output directory.

#### 3. Per-sample variant calling
`GATK HaplotypeCaller` will use BAM files as input for calling short variants (SNVs and INDELs) in GVCF mode. Default options are set as minimum base quality 20 and ploidy 1. This step will generate an intermediate file (GVCF) per-sample, which will be used for joint genotyping to produce a final multi-sample VCF file. A list of GVCF files generated will be parsed into tab-delimited text file (named as sample_map.txt) in the input directory. By default snpplet will proceed to the next step, all samples will be imported into the workspace for joint genotyping. You can also choose to stop snpplet at this step by specifying --stop true in snpplet.sh file. You can check the results of MultiQC (HTML reports) and Samtools (coverage.tsv) to select only samples that meet the criteria (e.g. minimum read quality, minimum mean depth, minimum genome coverage) before joint-genotyping. If this parameter is left as default setting, snpplet will proceed to step 4.

#### 4. Consolidating GVCFs for joint-genotyping
`GATK GenomicsDBImport` will create a GenomicsDB workspace and then import single-sample GVCFs into this workspace. GATK GenotypeGVCFs then will perform joint genotyping on GenomicsDB workspace to produce a multi-sample VCF. By default, `GenomicsDBImport` will read sample_map.txt to import all samples into the workspace. In case users have stopped the run after finishing step 3 for checking the results from steps 1 and 2 to exclude some problematic or low quality samples before proceeding to next steps. If you want to proceed to the next step without any change, use command ./sh snpplet.sh to do so. But if users want to include only some samples, users will need to provide a tab-delimited text file in input directory (formatted as shown below and file must be named sample_map_usr.txt)

			Sample_1 Sample_1.g.vcf.gz
			Sample_2 Sample_2.g.vcf.gz
			Sample_3 Sample_3.g.vcf.gz
                    .         .
                    .         .
                    Sample_n  Sample_n.g.vcf.gz

`GenomicsDBImport` then will automatically read sample_map_usr.txt file and import only samples that are listed into the workspace. After saving sample_map_usr.txt file to the input directory, execute command ./sh snpplet.sh to start joint genotyping.
										
#### 5. Variant filtering
By default, INDELS will be filtered with `GATK SelectVariants` (-select-type SNP). Remaining SNVs will be filtered by `GATK VariantsFiltration` with filter settings: QD < 2.0 || MQ < 40.0 (Filtered VCF file will be generated in vcf_joint directory in output directory (named as joint_filtered.vcf.gz)

#### 6. VCF to FASTA
Filtered VCF will be parsed to produce a SNV table using `GATK VariantsToTable`. SNV table will be parsed into FASTA format (file named aln.fasta) in aln directory in output directory. Position of all SNVs (in reference coordiates) will also be parsed into a separate file named as pos.txt in the same directory.



### Table 1 - Dependencies
Name | 
------ | 
Nextflow 20.07 (or later) | 
Java 8 | 
Fastqc | 
Multiqc | 
Trimmomatic | 
BWA | 
Samtools | 
Bcftools | 
Gatk 4.1.8.1 (or later) | 
Datamash | 

### Table 2 - Paths and parameters
Paths/Paramter | Paramter name | Default
------ | ------ | ------
`params.reads` | the location of the read fastq files | `/home/reads/*_{1,2}.fastq.gz`
`params.ref_genome` | the location of the read fastq files | `/home/ref/reference.fasta`
`params.results` | the location where the results will be stored |  `/home/results/`
`params.genome_name` | genome name | `NC000962.3`
`params.trimming_option` | parameter for trimming (Trimmomatic) | `SLIDINGWINDOW:4:30 MINLEN:70`
`params.mapping_option` | parameter for read mapping (BWA) | `-c 100 -M -T 50`
`params.haplotypecaller_option` | parameter for per-sample variant calling (HaplotypeCaller) | `-ploidy 1 -mbq 20`
`params.genomicsdbimport_option` | parameter for joint (GenomicsDBImport) | `test`
`params.genotypegvcfs_option` | parameter for joint (GenotypeGVCFs) | `test`
`params.variant_filter` | parameter for variant filtering |  `QD < 2.0` `MQ < 40.0`
`params.variant_filter_name` | parameter for variant filtering (VariantsFiltration) | `qd-mq`
`params.selectvariant_option` | parameter for variant filtering (SelectVariants) | `--exclude-filtered -select-type SNP`
`params.stop` | optional stop after step 3 (per-sample variant call) | `false`

## Quick installation of `snpplet`'s dependencies
#### 1. Install `homebrew` and `homebrew/science/bio` (bioinformatic fomulae)
<br>     1.1 Install `curl` and `git`:
          ``` 
          $ sudo apt-get update 
          ```
<br>     1.2 Install `homebrew`:
         ``` 
         $ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)" 
         ```
<br>     1.2.1 Set PATH:
         ``` 
         $ echo 'eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv)' >> /home/bb_18/.profile 
         ```
<br>     1.2.2 Evaluate `.profile`:
         ``` 
         $ eval $(/home/linuxbrew/.linuxbrew/bin/brew shellenv) 
         ```
<br>     1.3 Install `build-essential`:
          ``` 
          $ sudo apt-get install build-essential 
          ```
<br>     1.4 Install `gcc`:
          ``` 
          $ brew install gcc 
          ```
<br>     1.5 Tap `brewsci/bio`:
          ``` 
          $ brew tap brewsci/bio 
          ```
     
#### 2. Install dependencies for `snpplet` ( using `homebrew`) by using the following command:

     
     $ brew install nextflow \
     > gatk \
     > fastqc \
     > samtools \
     > bwa \
     > datamash \
     > bcftools \
     > trimmomatic \
     > samtools
     

   <br>  2.1 Launch `nexflow` to install its dependencies by running `$ nextflow`
     
#### 3. Install another dependency for `snpplet` using `pip`:
        
        $ pip3 install multiqc
        





# Schematic outline
![Image of Yaktocat](https://github.com/CENMIG/snpplet/blob/master/figures/tuber.png)

