# tuber
## Description
Whole genome sequencing data are used for various purposes ///Variant detection, phylogenetic analysis, …///. This pipeline is an implementation of /a series of/ tools for pre-processing raw data, quality control checking, mapping (against the reference genome), (short) variant calling and filtering to produce a ready to use vcf/msa file. The GATK Best Practices implemented in this pipeline are used for variant calling in GVCF mode, which allows users to work on a large set of samples. 

tuber is composed of two parts, ////preprocessing+per-sample variant calling  ///joint-genotyping+VCF-to-MSA.

Part I (tuber_1.nf)
preprocessing-mapping-sorting-per sample variant calling
Results of this part will help users to QC check and select a set of high quality sample for joint-genotyping.
Users will need to create a sample map file 

Part II (tuber_2.nf)
flexibility to 

This pipeline applies workflow 
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

|  | Part I | Part II |
| ------ | ------ | ------ |
| Input | fastq, fasta| GVCFs, fasta, sample map
| Output | GVCFs, HTML, Table | VCF, MSA (.fasta)

Raw data in fastq format, file name should match the following scheme:  _{1,2}.gz. Reference genome sequence in fasta format. Tuber will produce outputs inside the directory specified by users. Final output will be filtered VCF file (.vcf) and multiple sequence alignment (.fasta).

# Part I

 ```sh
     $ nextflow run >>>>>1.nf
```

#### This script will perform the following:
#### 1. Preprocessing raw reads
Trimmomatic will be used for Illumina adapters clipping and trimming low quality (sequences in the) reads. Short reads will be dropped. By default it is set as: 
- PE -phred33, ILLUMINACLIP:TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:4:30, MINLEN:70

Fastqc is used to QC check reads before and after processing with Trimmomatic. Then MultiQC will parse summary stats from fastqc results of each sample to produce HTML report and summary tables of multiple sample.

#### 2. Mapping
Processed reads will be aligned to the reference genome sequence with BWA. Default is mem -c 100 -M -T 50
 SAM files will be converted into BAM files by samtools and duplicate reads will be marked by Picard.
Average depth and genome coverage are calculated using Samtools.

#### 4. Per-sample variant calling
HaplotypeCaller will use BAM files as input to call short variants (SNVs and INDELs) in GVCF mode. Minimum base quality is set to 20. This step will generate an intermediate file (GVCF) per-sample, which will be used for joint genotyping to produce a final multi-sample VCF file.

#### Generating a sample map file (by users)
Combine results processes 1b and 2c, samples that meet the criteria will be listed in a sample map file. This file will be used in 4a (consolidating GVCFs). Default is set to 1) GC content (passed) 2) meandepth>10X 3) coverage>90%

# Part II 
Users will need to create a sample map file (file name must be sample_map.txt), save this to ____ directory.
Then use this command to run part 2 of pipeline.

  ```sh
     $ nextflow run >>>>>2.nf
```
#### 1. Consolidating GVCFs and joint-genotyping
GATK GenomicsDBImport will create a GenomicsDB workspace and then import single-sample GVCFs (listed in a sample map file) into this workspace. GATK GenotypeGVCFs then will perform joint genotyping on GenomicsDB workspace to produce a final VCF.
#### 2. Hard filtering and annotation
INDELS will be filtered with GATK SelectVariants. Remaining SNVs will be filtered by GATK VariantsFiltration with default filter settings 1) minimum QD of 2 and 2) minimum MQ of 50.

#### 3. VCF to FASTA
Filtered VCF will be pasred to produce SNV table using GATK VariantsToTable.
SNV table will be parsed into FASTA format.

# Installing and usage
1. Install dependencies listed below
2. Download nextflow pipeline scripts to (user specified directory)
     ```sh
     $ curl https://github.com/b600a/<<<<<<<.nf
     $ curl https://github.com/b600a/>>>>>>>.nf
     ```
3. Users will need to specify 1) the location of raw read (FASTQ) files, 2) the location of reference genome sequence (FASTA) file and 3) the location of output. Users will also need to specify some parameters to make them more suitable for the dataset. By using a text editor, look sections listed below in <<<<<<<.nf file to define proper parameters.

    | Section | Parameter |
    | ------ | ------ |
    | User Defines LOCATIONS | The location of raw read (FASTQ) files, reference genome sequence (FASTA) file, output directory|
    | User Defines TRIMMING OPTIONS | phred quality score, adapter clipping, SLIDINGWINDOW:4:30, minimum read lenght |
    | User Defines SAMPLE FILTERING OPTIONS | minimum mean depth, minimum genome coverage |
    | User Defines HARD FILTERING OPTIONS | QD, MQ |


4. Launch the pipeline execution with the following command:
    ```sh
    $ nextflow run >>>>>>>>>>>>.nf
    ```
# Schematic outline
insert picture here

# Dependencies

Nextflow
	https://www.nextflow.io/docs/latest/getstarted.html#installation
Fastqc 
	https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Multiqc
	https://multiqc.info/, https://github.com/ewels/MultiQC
Trimmomatic
	http://www.usadellab.org/cms/?page=trimmomatic
BWA
	https://github.com/lh3/bwa
Samtools

Bcftools
	
GATK4
https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4




[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>
