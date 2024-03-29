#!/bin/bash
# This script extract the reads from the given BAM file (mapped to GRCh38 reference) using input BED regions 
# and then remap the reads to the GRCh38 reference sequences that contains the decoy sequences.  

THREAD=1 #set number of threads Maybe later over parameter!? 

usage="sh ShortReadRemap.sh [-h] [-b -d -r -o] -- script to extract the reads from the given bed regions and remapping them to reference with decoy sequences

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
    echo "sh ShortReadRemap.sh [-h] -b <bam/cram> -d <bed> -r <reference> -o <out_prefix>"
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

# Extract the BAM corresponds to the input bed regions and excludes reads unmappd, mate unmapped, not primary alignment, and supllementary alignment (-F 2316)
	echo -e "\tExtract reads for remapping"
        if [[ $bam =~ \.cram$ ]]
        then
	    samtools view -F 2316 -hb --reference ${old_ref} -o ${out}_extracted_reads.bam ${bam} `perl -ane '{print "$F[0]:$F[1]-$F[2] "}' $bed `
        else
            samtools view -F 2316 -hb -o ${out}_extracted_reads.bam ${bam} `perl -ane '{print "$F[0]:$F[1]-$F[2] "}' $bed `
        fi
# Read names unique in file and extract BAM corresponds to paired-end reads
	samtools view -hb -F 2316 ${out}_extracted_reads.bam | samtools sort -@ ${THREAD} -n -  > ${out}_extracted_reads_sorted.bam 
	samtools view ${out}_extracted_reads_sorted.bam  | cut -f1 | uniq -c | awk '$1!=2 {print $2}' > ${out}_nonpairs_rnames.txt
	samtools view -Sh ${out}_extracted_reads_sorted.bam | fgrep -vf ${out}_nonpairs_rnames.txt | samtools view -hb - > ${out}_pairs_only.bam 
        
	samtools sort -@ ${THREAD} -n -o ${out}_original_sorted_by_read_names.bam ${out}_pairs_only.bam
	samtools sort -@ ${THREAD} -o ${out}_original_sorted.bam ${out}_pairs_only.bam
	samtools index ${out}_original_sorted.bam
	
	echo -e "\tCovert reads to fastq"
# Extract the paired-end reads
	samtools view -H ${out}_extracted_reads.bam | grep "^@RG"| sed 's/	/\\t/g'  > ${out}_bug_regions.RG.txt
        samtools fastq -@ ${THREAD} ${out}_original_sorted_by_read_names.bam -1 ${out}_extract_1.fastq -2 ${out}_extract_2.fastq -0 /dev/null -s /dev/null -n
	
	echo -e "\tStart remapping"
# Map the extracted reads to the reference 
	bwa mem -t ${THREAD} -R $(head -n 1 ${out}_bug_regions.RG.txt) ${ref} ${out}_extract_1.fastq ${out}_extract_2.fastq | samtools view -hb - > ${out}_remapped.bam
	samtools sort -@ ${THREAD} -o ${out}_remapped_sorted.bam ${out}_remapped.bam 
	samtools index ${out}_remapped_sorted.bam


# Remove intermediate/temp files
	rm ${out}_remapped.bam
	rm rm ${out}_extracted_reads.bam
        rm ${out}_extracted_reads_sorted.bam
	rm ${out}_original_sorted_by_read_names.bam
	rm ${out}_pairs_only.bam
        rm ${out}_nonpairs_rnames.txt
        rm ${out}_bug_regions.RG.txt
        rm ${out}_extract_2.fastq
        rm ${out}_extract_1.fastq
	echo "Scucessfully finsihed"
	#exit 0
else

	echo "Error one or more files do not exist. Check your paths"
	exit 1
fi

echo "<out_prefix>_original_sorted.bam <- Sorted BAM files corresponds to the extracted paired-end reads (Original mapping)"
echo "<out_prefix>_remapped_sorted.bam <- Sorted bam file after re-mapping extracted reads to GRCh38 masked v2 reference (mapping to modified reference)"
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "Script Execution Time: $execution_time"
