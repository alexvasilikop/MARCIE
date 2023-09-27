#!/bin/bash

#exit immediately if a command returns non-zero status
set -eo pipefail

#Check that required software is in path
for software in bwa gatk sambamba picard samtools; do

    if ! command -v $software &> /dev/null; then
        echo -e "### $software could not be found ###\n"
        return 127
    else
        echo -e "\n### $software is on path ###\n"

    fi
done

#Parse input variables from command line
#genome
genome=$1
echo -e "\n### Genome file provided: $genome"

#reads
READS_DIR=$2
echo -e "\n### Directory with NGS reads (Illumina) for the different samples: $READS_DIR"

#species_name
species_name=$3
echo -e "### Species name for RG string: $species_name"

#date
current_date=$(date +%F)
echo -e "### Date for RG string: $current_date"

#Threads for mapping and processing of BAM files
threads=$4
echo -e "### Number of threads to use: $threads"

#Minimum mapping quality for the reads to keep (lower quality mappings are filtered out)
MQ=$5
echo -e "### Minimum mapping quality to filter BAM (MQ): $MQ"

echo -e "\n### Sample IDs (will be used for RG string in the BAM):"
ls -1 $READS_DIR/*fastq.gz | xargs -n1 basename | cut -f1 -d"." | sort -k1 -n | uniq > list_of_SeqIDs.txt
cat list_of_SeqIDs.txt

############################################################################
#Step 1: Read mapping with BWA and processing of BAM files for variant calling
echo -e "\n### Step 1: Read mapping with BWA and processing of BAM files for variant calling... ##\n"
###################################################################################################

#perform mapping of reads and filtering of bam files
cat list_of_SeqIDs.txt| while read SeqID;do

    mappingID=$(sed -e "s/.fasta//" <<< $genome).$SeqID
	
    #check if mapping files with marked duplicates already exist
	if [[ -f $mappingID.sorted.MQ$MQ.mdup.bam ]];then
        echo -e "### Mapping files for $SeqID exist... ###\n"
        continue

    else
        ###Perform mapping and processing of BAM files for variant calling)
        #checks for reads files
        echo -e "\n### Working on sample $SeqID###\n"
        FILE1="$SeqID"".R1.fastq.gz"
        echo -e "\n###fastq R1: $FILE1###"
        #check if file1 exists
        [[ -f "$READS_DIR/$FILE1" ]] && echo -e "### $FILE1 exists in reads directory."

        FILE2="$SeqID"".R2.fastq.gz"
        echo -e "\n###fastq R2: $FILE2###"
        #check if file 2 exists
        [[ -f "$READS_DIR/$FILE2" ]] && echo -e "### $FILE2 exists in reads directory.\n"

        #Create RGstring
        RGstring="@RG\tID:$SeqID-illumina_reads\tDT:$current_date\tLB:$SeqID-illumina_reads\tPL:ILLUMINA\tSM:$species_name""_ref_$SeqID-illumina_reads"
        echo -e "\n###RGstring:$RGstring"
        
        #indexing genome for BWA
        echo -e "\n### Indexing genome reference: ###\n\t$genome\n"
        bwa index $genome

        #Mapping reads and create mapping files (option -M for picard compatibility)
        echo -e "\n### Mapping Reads onto genome reference\ngenome:\t$genome\nfq1:\t$READS_DIR/$FILE1\nfq2:\t$READS_DIR/$FILE2\n"
        bwa mem -t $threads -M -R $RGstring $genome $READS_DIR/$FILE1 $READS_DIR/$FILE2 -o $mappingID.sam

        #Post-treating mapping files (SAM->BAM, sort, filter, mark duplicates)
        echo -e "\n### Post-treating mapping files\nmappingID:\t$mappingID\nSAM->BAM and sorting...\n"
        sambamba view -t $threads -f bam -h -S -o $mappingID.bam $mappingID.sam

        echo -e "\n### Deleting SAM file ...###\n"
        rm -f $mappingID.sam
	
	echo -e "\n### Sorting BAM file ...###\n"
        sambamba sort -t $threads --tmpdir=./ $mappingID.bam -o $mappingID.sorted.bam
        echo -e "\n### Deleting unsorted BAM file ...###\n"
        rm -f $mappingID.bam

        #filter reads with mapping quality <10
        echo -e "\n### Filtering reads with mapping quality < 10\nmappingID:\t$mappingID\nMinimal Mapping Quality:\t$MQ"
        sambamba view -t $threads -f bam -h -F "mapping_quality >= $MQ" $mappingID.sorted.bam -o $mappingID.sorted.MQ$MQ.bam

        #mark duplicates (necessary for variant calling)
        echo -e "\n### Mark duplicates with Picard\nmappingID:\t$mappingID\n"
        time picard MarkDuplicates -I $mappingID.sorted.MQ$MQ.bam -O $mappingID.sorted.MQ$MQ.mdup.bam -M marked_dup_metrics.$mappingID.sorted.MQ$MQ.txt

        #remove temporary files
        echo -e "\n### Remove temporary files (Deleting sorted filtered BAM file before quality-based filtering step) ...###\n"
        rm -f $mappingID.sorted.MQ$MQ.bam
        rm -f $mappingID.sorted.MQ$MQ.bam.bai
	rm -f $mappingID.sorted.bam
 	rm -f $mappingID.sorted.bam.bai

        #index final bam file
        echo -e "\n### Index mapping file (BAM file) \nmappingID:\t$mappingID:\n"
        sambamba index -t $threads $mappingID.sorted.MQ$MQ.mdup.bam
    fi
done


###############################
#Variant calling (step 2)
echo -e "\n### Step 2: Short germline variant calling with HaplotypeCaller... ##\n"
#####################################################################

cat list_of_SeqIDs.txt| while read SeqID;do

    mappingID=$(basename $genome .fasta).$SeqID
    #input bam for variant calling
    input=$mappingID.sorted.MQ$MQ.mdup.bam

    if [[ -f $input.g.vcf ]];then
        echo -e "### GVCF of sample $SeqID ($input) exists... ###\n"
        continue
    else
        echo -e "\n### Input BAM $mappingID.sorted.MQ$MQ.mdup.bam...\n"

        #Genome dictionary and indexing
        echo -e "\n### Create genome dictionary with picard...\ngenome:\t$genome\n"
        ref_genome_dict=$(basename $genome .fasta).dict

        if [ ! -f "$ref_genome_dict" ];then
            time picard CreateSequenceDictionary -R $genome -O $ref_genome_dict
        else
            echo -e "### Genome dictionary found in current directory... ###\n"
        fi

        #Genome indexing with samtools
        echo -e "\n### Index genome with SAMTOOLS...\ngenome:\t$genome\nbam file:\t$input"
        samtools faidx $genome

        #Call short germline variants (SNPs and indels) using HaplotypeCaller
        echo -e "\n### SNP-Calling with GATK \ngenome:\t$genome\nbam file:\t$input"
        time gatk HaplotypeCaller -R $genome -I $input -ERC GVCF -O $input.g.vcf
    fi
done


############################################################
#STEP 3: COMBINE GVCFS 
echo -e "\n### Step 3: Combine GVCFs with \"gatk CombineGVCFs\"... ###\n"
#############################################################################

list_GVCFs=$species_name".MQ"$MQ".mdup.bam.g.vcf.list"
combined_gvcfs=$species_name".combined.g.vcf"
if [[ -f $combined_gvcfs ]]; then
    echo -e "\n### File with combined GVCFs already exists... ###\n"
else
    echo -e "\n### List of GVCFs: $list_GVCFs ###\n"
    #WARNING: IF the file has a different ending than *list the CombineGVCFs command gives an error
    ls -1 *g.vcf > $list_GVCFs
fi

#Combine gvcfs
time gatk CombineGVCFs -R $genome -V $list_GVCFs -O $combined_gvcfs


###############################################################################
## STEP 4: GENOTYPING COMBINED GVCFS - GATK  ##################
echo -e "\n### Step 4: Genotyping combined GVCFs with \"gatk GenotypeGVCFs\"... ###\n"
###########################################################################################

#make output directory for placing the final combined file with all genotypes from diffrent samples
outdir=joint_vcf_raw.$genome
raw_vcf="$species_name.raw.vcf"
mkdir -p $outdir

if [[ -f "$outdir/$raw_vcf" ]]; then
    echo -e "\n### Combined file with genotyped GVCFs already exists... ###\n"
else
    time gatk GenotypeGVCFs -R $genome -V $species_name.combined.g.vcf -O $outdir/$raw_vcf
fi

#Pipeline finished
echo -e "\n### Variant calling pipeline finished succesfully... ###\n"
