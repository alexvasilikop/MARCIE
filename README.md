# MARCIE (Mapping, vARiant Calling and genotyping pIpEline)

## Introduction
This is a simple pipeline written in BASH for performing reads mapping (Illumina paired-end reads) on a reference genome 
and subsequent filtering of BAM alignment files, short germline variant calling and combined genotyping of several 
samples (using GATK).

Specifically the steps are the following:
1. Check for software dependencies
2. Run mapping of all libraries using the BWA mem algorithm in paired-end mode
3. Converting, sorting and filtering of BAM files according to provided threshold for mapping quality (MQ) - each library is processed searately after each mapping is finished
4. Marking duplicates with Picard
5. Short germline variant calling with HaplotypeCaller (GATK toolkit)
6. Combining the GVCF files (CombineGVCFs)
7. Genotyping the combined GVCF files (GenotypeGVCFs)

## Requirements
1. BWA mapper
2. GATK toolkit (for variant calling and subsequent genotyping)
3. Sambamba tool for processing and filtering BAM files (sorting, mapping quality filtering)
4. Picard (for marking duplicates)
5. Samtools (for indexing the genome fasta files)

The script requires that all above software dependencies are on the path and performes an initial check for dependencies.
The easiest way to run the pipeline is to generate a conda environment, install all dependencies there and 
run the pipeline while the environment is activated.

```
#create and activate env
conda create -n mapping_varcall_genotype
conda activate mapping_varcall_genotype
#install requirements
conda install -c bioconda sambamba
conda install -c bioconda samtools
conda install -c bioconda picard
conda install -c bioconda bwa
conda install -c bioconda gatk
```

## Formatting of input files
The script takes as input a list of paired-end Illumina files representing different samples or libraries (e.g., see below `sample_1`, `sample_2`) and a genome assembly file in fasta 
format. A directory path with all the reads' fastq files has to be provided as a command line argument to the bash script. Make sure the samples/libraries have the following naming format 
in your reads' directory (replace `sample_1`, `sample_2` with your sample codes). The pipeline assumes gzipped fastq files of the Illumina reads. If another format or naming is used 
the pipeline will fail. 

Example 1:
```
sample_1.R1.fastq.gz
sample_1.R2.fastq.gz
sample_2.R1.fastq.gz
sample_2.R2.fastq.gz
```

Example 2:
```
codeA.R1.fastq.gz
codeA.R2.fastq.gz
codeB.R1.fastq.gz
codeB.R2.fastq.gz
```
## Usage
Provide the following inputs to the ´marcie.sh´ script with the following order (the order is important but the variable names might differ. You could also provide the paths to these 
files as command line arguments to the script without assigning them to variables beforehand):

Example:
```
#set variables
genome=example_genome.fasta
READS_DIR=/home/user/reads
SPECIES=species_sp
THREADS=10
#threshold for filtering out low quality mappings
MQ=10

#run pipeline
./marcie.sh $genome $READS_DIR $SPECIES $THREADS $MQ
```

The species name variable is important as it will be the prefix of all your output files.

