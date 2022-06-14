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
        s) sequence=${OPTARG};;
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
        if [[ $bam =~ \.cram$ ]]
        then
            samtools view -F 2316 -hb --reference ${old_ref} -o ${out}_extracted.bam ${bam} `perl -ane '{print "$F[0]:$F[1]-$F[2] "}' ${bed} `
        else
            samtools view -F 2316 -hb -o ${out}_extracted.bam ${bam} `perl -ane '{print "$F[0]:$F[1]-$F[2] "}' ${bed} `
        fi
	
	echo -e "\tCovert reads to fastq"
# extract the paired-end reads
	samtools bam2fq ${out}_extracted.bam > ${out}_extracted.fastq
  #samtools fastq -0 /dev/null -s /dev/null -n -@ ${THREAD} ${out}_extracted.bam > ${out}_extracted.fastq
	
	echo -e "\tStart remapping"
# map the extracted reads to the reference 
  if [[ $sequence == "ONT" ]]
    minimap2 -ax map-ont ${ref} ${out}_extracted.fastq > ${out}_remapped.sam
    samtools view -bS ${out}_remapped.sam > ${out}_remapped.bam
    samtools sort ${out}_remapped.bam > ${out}_remapped.sorted.sam
    rm ${out}_extracted.fastq
    rm ${out}_remapped.sam
    rm ${out}_remapped.bam
  elif [[ $sequence == "PB" ]]
	  minimap2 -ax map-pb ${ref} ${out}_extracted.fastq > ${out}_remapped.sam
    samtools view -bS ${out}_remapped.sam > ${out}_remapped.bam
    samtools sort ${out}_remapped.bam > ${out}_remapped.sorted.sam
    rm ${out}_extracted.fastq
    rm ${out}_remapped.sam
    rm ${out}_remapped.bam
	echo "Scucessfully finsihed"
	exit 0
else

	echo "Error one or more files do not exist. Check your paths"
	exit 1
fi

echo "<out_prefix>_extracted.bam <- Sorted BAM files corresponds to the extracted paired-end reads (original mapping)"
echo "<out_prefix>_remapped.sorted.bam <- Sorted bam file after minimap step i.e. mapping extracted reads to GRCh38 masked v2 reference (final output)"
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "Script Execution Time: $execution_time"
