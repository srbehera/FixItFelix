#!/bin/bash
# This script extract the reads from the given BAM file (mapped to GRCh38 reference) using input BED regions 
# and then remap the reads to the GRCh38 reference sequences that contains the decoy sequences.  

THREAD=1 #set number of threads Maybe later over parameter!? 

usage="sh script.sh [-h] [-b -d -r -o] -- script to extract the reads from the given bed regions and remapping them to reference with decoy sequences

where:
    -h  show this help text
    -b  path to the BAM/CRAM file
    -f  Old reference (if CRAM file is used)
    -d  path to the BED file
    -r  Reference file with decoy sequences
    -o  output prefix
    -t  number of threads"

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` $usage"
  exit 0
fi

if [ $# -eq 0 ]; then
    echo "Error: No arguments provided"
    echo "sh script.sh [-h] -b <bam/cram> -d <bed> -r <reference> -o <out_prefix>"
    exit 1
fi

while getopts b:d:r:o:f:t: flag
do
    case "${flag}" in
        b) bam=${OPTARG};;
        d) bed=${OPTARG};;
        r) ref=${OPTARG};;
        o) out=${OPTARG};;
        f) old_ref=${OPTARG};;
        t) THREAD=${OPTARG};;
    esac
done

if [[ $bam =~ \.cram$ &&  -z "$old_ref" ]]; then
    echo "Error: For CRAM files, use -f parameter with old reference"
    exit 1
fi

echo "BAM/CRAM file: $bam";
echo "Bed file: $bed";
echo "Reference: $ref";
echo "Output: $out";
echo "Old Ref: $old_ref";

##cd ${PWD}
# cd ${PBS_O_WORKDIR}
start=$(date +%s.%N)

if [[ -e $ref && -e $bam && -e $bed ]] ; then 
##start=$(date +%s.%N)

# extract the bam corresponds to the bed region
	echo -e "\tExtract reads for remapping"
##samtools view -hb -L ${bed} -o ${out}/bug_regions.bam ${bam}
        samtools fastq -@ ${THREAD} ${bam} -1 ${out}.R1.fastq.gz -2 ${out}.R2.fastq.gz -0 /dev/null -s /dev/null -n
	STAR --runMode alignReads --readFilesCommand zcat --genomeDir <index_file_modified_GRCh38> --outFileNamePrefix ${out}_remapped --readFilesIn ${out}.R1.fastq.gz ${out}.R2.fastq.gz
	echo "Scucessfully finsihed"
	exit 0
else

	echo "Error one or more files do not exist. Check your paths"
	exit 1
fi

#echo "<out_prefix>_pairs_only.bam <- Sorted BAM files corresponds to the extracted paired-end reads (original mapping)"
echo "<out_prefix>_remapped.bam <- STAR remapped BAM file"
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "Script Execution Time: $execution_time"
