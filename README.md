# FixItFelix: Fixing GRCh38 bugs
Investigation of collapsed and duplicated regions in the GRCh38 (known as `hg38`) reference sequence

This is a script and resources to realign short paired-end reads (also long reads) across medical relevant but challenging regions. The bed file holds the regional coordinates from GRCh38 regions that are either falsely duplicated or collapsed. As such these errors impact ~4.1Mbp of sequence that can be improved using this remapping step.

The tool has three different scripts for short read (illumina), long read (ONT or PacBio) and RNA-Seq datasets. For long reads, we are using [`minimap2`](https://github.com/lh3/minimap2) aligner, so the parameters need to be changed depending on sequencing technologies. For example, `-ax map-ont` for ONT reads and `-ax hifi` for HiFi reads. For RNA-Seq data, we extract all the sequences from the BAM file instead of regions of given BED file. We used [`STAR`](https://github.com/alexdobin/STAR) aligner for re-aligning the extracted sequences to the modified reference genome. Please change the command if [`HISAT2`](https://github.com/DaehwanKimLab/hisat2) or other splice-aware aligners are used. Also, please generate the index file for the reference and put the path in the command.

## Running the script
To run the script you need a modified [GRCh38 reference from here](https://bcm.box.com/s/xi95ahgzrw86pvogm7sdwl0ppn49i5dn). We also provided [2nd version of modified GRCh38 reference](https://bcm.box.com/s/qz5h36ry4cg9j15mwolzcpf1z7ha823e) that excludes decoys related to GPRIN2, DUSP22 and FANCD2 genes. Please use the bed file in this repository (add +/- 5kbp flanking regions). The script assumes a global installation of `samtools`, `bwa`, `minimap2`, `STAR` and all required softwares. Furthermore, please have the reference index for `bwa` (see [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml))  at the same location as the reference! 

### For WGS or WES dataset
    [Short/Long]ReadRemap.sh] -b <BAM> -r <Modified_GRCh38_Ref> -d <bed_regions> -o <out_prefix>
    
 If CRAM file is used instead of BAM, then it is required to use the [original GRCh38 reference](https://bcm.box.com/s/ym4x3z61ib4okbguy7zn8lre8uo6mcxz) in the script with `-f` parameter.
 
    [Short/Long]ReadRemap.sh -b <BAM> -r <Modified_GRCh38_Ref> -d <bed_regions> -f <Original_GRCh38_Ref> -o <out_prefix>
    
 To run with multiple threads, use `-t <no_threads>`
 
 ### For RNA-Seq datase
 
    RNASeqRemap.sh -b <BAM> -r <Modified_GRCh38_Ref> -o <out_prefix>
    
## Output

The scripts produce following two output file 

    1. <out_prefix>_original_sorted.bam <- Sorted BAM files corresponds to the extracted paired-end reads (Original mapping)"
    2. <out_prefix>_remapped_sorted.bam <- Sorted bam file after BWA MEM step i.e. re-mapping of extracted reads to modified GRCh38 reference (final output)" 
    
The above two BAM files are corresponding to the BED regions as it realigns only the regions that are impacted due to falsely duplicated or collapsed errors. For getting a BAM file for the whole genome, it needs to be merged with the mapping of unaffected regions. The current version of the script does not have a merging step at the end.
