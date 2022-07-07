# FixItFelix: Fixing GRCh38 bugs
Investigation of collapsed and duplicate regions in the HG38 reference sequence

This is a script and resources to realigning short paired-end reads across medical relevant but challenging regions. The bed file holds the regional coordinates from GRCh38 regions that are either wrongly duplicated or collapsed. As such these errors impact ~4.1Mbp of sequence that can be improved using this remapping step.


## Running the script
To run the script you need a modified [GRCH38 reference from here](https://bcm.box.com/s/xi95ahgzrw86pvogm7sdwl0ppn49i5dn). Furthermore use the bed file in this repository. The script assumes a global installation of `samtools` and `bwa`. Furthermore have the reference index for `bwa` (see [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml))  at the same location as the reference! 

    sh extract_remap.sh -b <BAM> -r <Reference_with_decoy> -d <bed_regions> -o <out_prefix>
    
 If CRAM file is used instead of BAM, then it is required to use the [old reference](https://bcm.box.com/s/ym4x3z61ib4okbguy7zn8lre8uo6mcxz) in the script with `-f` parameter.
 
    sh extract_remap.sh -b <CRAM> -r <Reference_with_decoy> -d <bed_regions> -f <old_reference> -o <out_prefix>
    
 To run with multiple threads, use `-t <no_threads>`
    
## Output

    <out_prefix>_original.sorted.bam <- Sorted BAM files corresponds to the extracted paired-end reads (original mapping)"
    <out_prefix>_remapped.sorted.bam <- Sorted bam file after BWA MEM step i.e. re-mapping of extracted reads to GRCh38 masked v2 reference (final output)"
