# tuber
# Description
Whole genome sequencing data are used for various purposes including molecular typing, classification and phylogenetic analysis./ We present ‘tuber’, a Nextflow pipeline for variant calling from WGS data. This pipeline is an implementation of tools for pre-processing raw data and quality control, mapping, (short) variant calling and filtering to produce a ready to use vcf and multiple sequence alignment files. tuber includes workflow of GATK Best Practices for germline short variant discovery, [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). VQSR is not included here and hard filtering is used instead. In addition to the GATK Best Practices, tuber includes steps to parse the SNV table into a multiple sequence alignment in fasta format.

# Inputs
- raw data (.fastq); file name should/must match the following scheme: _{1,2}.fastq.gz 
- reference genome sequence (.fasta)
- (Optional) adapter sequence (.fasta); learn more at [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- (Optional) list of sample id to be included for joint-genotyping (must be named sample_map_usr.txt)

# Outputs 
tuber will produce outputs inside the output directory specified by users.
-HTML report of FastQC 
-tables, many tables
- A final VCF
- multiple sequence alignment in fasta format

# Installation and usage
1. Install dependencies listed below. (tuber will run on Linux or MacOSX; )
2. Download nextflow pipeline scripts to (user specified directory)
     ```sh
     $ git clone https://github.com/b600a/tuber.git
     $ cd tuber
     ```
3. Users will need to specify paths to the input files, reference genome sequence and path for outputs by using any text editor to edit tuber.sh file (--ref_genome /////, /////). Users may also need to specify some parameters to make them more suitable for the dataset, these are listed below. 

    | Section | Parameter | Default |
    | ------ | ------ | ------ |
    | paths & inputs | `baseDir` <br>`params.ref_genome` <br>`params.reads` <br>`params.results` <br>`params.genome_name` | 
    | parameters for trimming | `params.trimming_option` |
    | parameter for read mapping | `params.mapping_option` |
    | parameters for variant calling | `params.haplotypecaller_option` <br>`params.genomicsdbimport_option` <br>`params.genotypegvcfs_option` |
    | parameters for variant filtering | `params.variant_filter` <br>`params.variant_filter_name`  <br>`params.selectvariant_option` |
    | optional stop | params.stop | false |

4. Launch the pipeline execution with the following command:
    ```sh
    $ ./run.sh
    ```
5. tuber will consume huge amount of space. Once the run has finised, you can clean up workspace by removing `work/`. 

#### This script will perform the following:
#### 1. Preprocessing raw reads
Trimmomatic will be used for Illumina adapters clipping and trimming low quality (sequences in the) reads. Reads that are shorter than the defined length will be dropped. By default it is set as: `-PE -phred33, ILLUMINACLIP:TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:4:30, MINLEN:70`

Fastqc will be used forQC checking sequencing reads before and after processing with Trimmomatic. Then MultiQC will parse summary stats from fastqc results of each sample to produce summary tables (multiqc_fasqc.txt and multiqc_general_stats.txt and HTML reports in `multiqc_fastqc/` in output directory.

#### 2. Mapping, indexing, sorting, calculating depth and coverage
Processed reads will be aligned to the reference genome sequence with BWA MEM. Default is set as: `-c 100 -M -T 50`. SAM files will be converted into BAM files by samtools and duplicate reads will be marked by Picard. Average depth and genome coverage of each sample are calculated from BAM using Samtools, results will be parsed into coverage.tsv file in `bam/` in output directory.

#### 3. Per-sample variant calling
GATK HaplotypeCaller will use BAM files as input for calling short variants (SNVs and INDELs) in GVCF mode. Default options are set as Minimum base quality 20 and ploidy 1. This step will generate an intermediate file (GVCF) per-sample, which will be used for joint genotyping to produce a final multi-sample VCF file. A list of GVCF files generated will be parsed into tab-delimited text file (named as sample_map.txt) in the input directory. By default tuber will proceed to the next step, all samples will be imported into the workspace for joint genotyping. Users can also choose to stop tuber at this step by specifying `--stop true` in tuber.sh file. Users can check the results of MultiQC (HTML reports) and Samtools (coverage.tsv) to select only samples that meet the criteria (e.g. minimum read quality, minimum mean depth, minimum genome coverage) before joint-genotyping. If this parameter is left as default setting, tuber will proceed to step 4.

#### 4. Consolidating GVCFs and joint-genotyping
GATK GenomicsDBImport will create a GenomicsDB workspace and then import single-sample GVCFs into this workspace. GATK GenotypeGVCFs then will perform joint genotyping on GenomicsDB workspace to produce a multi-sample VCF. By default, GenomicsDBImport will read sample_map.txt to import all samples into the workspace. In case users have stopped the run after finishing step 3 for checking the results from steps 1 and 2 to exclude some /problematic, low quality/ samples before proceeding to next steps. If users want to proceed to the next step without any change, use command `./sh tuber.sh` to do so. But if users want to include only some samples, users will need to provide a tab-delimited text file in input directory (formatted as shown below and file must be named sample_map_usr.txt)

Sample_1	Sample_1.g.vcf.gz
Sample_2	Sample_2.g.vcf.gz
Sample_3	Sample_2.g.vcf.gz

GenomicsDBImport then will automatically read sample_map_usr.txt file and import only samples that are listed to the workspace. After saving sample_map_usr.txt file to the input directory, execute command ./sh tuber.sh to start joint genotyping.

#### 5. Variant filtering
By default, INDELS will be filtered with GATK SelectVariants. Remaining SNVs will be filtered by GATK VariantsFiltration with filter settings: minimum QD of 2 and minimum MQ of 40. Filtered VCF file will be generated in vcf_joint directory in output directory (named joint_filtered.vcf.gz)

#### 6. VCF to FASTA
Filtered VCF will be parsed to produce a SNV table using GATK VariantsToTable. 
SNV table will be parsed into FASTA format (file named aln.fasta) in aln directory in output directory.
Position of SNVs (in reference coordiates) will also be parsed into file named pos.txt in the same directory.



# Schematic outline
insert picture here

# Dependencies

* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
* [Fastqc ](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )
* [Multiqc](https://multiqc.info/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic )
* [BWA](https://github.com/lh3/bwa )
* [Samtools](http://www.htslib.org)	
* [Bcftools](http://www.htslib.org)
* [GATK4](https://github.com/broadinstitute/gatk/releases)
