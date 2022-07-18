# FixItFelix: Fixing GRCh38 bugs
Investigation of collapsed and duplicated regions in the GRCh38 (known as `hg38`) reference sequence

This is a script and resources to realigning short paired-end reads across medical relevant but challenging regions. The bed file holds the regional coordinates from GRCh38 regions that are either wrongly duplicated or collapsed. As such these errors impact ~4.1Mbp of sequence that can be improved using this remapping step.

The tool has three different scripts for short read (illumina), long read (ONT or PacBio) and RNA-Seq datasets. For long reads, we are using [`minimap2`](https://github.com/lh3/minimap2) aligner, so the parameters need to be changed depending on sequencing technologies. For example, `-ax map-ont` for ONT reads and `-ax hifi` for HiFi reads. For RNA-Seq data, we extract all the sequences from the BAM file instead of regions of given BED file. We used [`STAR`](https://github.com/alexdobin/STAR) aligner for re-aligning the extracted sequences to the modified reference genome. Please change the command if [`HISAT2`](https://github.com/DaehwanKimLab/hisat2) or other splice-aware aligners are used. Also, please generate the index file for the reference and put the path in the command.

## Running the script
To run the script you need a modified [GRCh38 reference from here](https://bcm.box.com/s/xi95ahgzrw86pvogm7sdwl0ppn49i5dn). We also provided [2nd version of modified GRCh38 reference](https://bcm.box.com/s/qz5h36ry4cg9j15mwolzcpf1z7ha823e) that excludes decoys related to GPRIN2, DUSP22 and FANCD2 genes. Furthermore use the bed file in this repository (add +/- 5kbp flanking regions). The script assumes a global installation of `samtools`, `bwa`, `minimap2`, `STAR' and all required softwares. Furthermore have the reference index for `bwa` (see [bwa index](http://bio-bwa.sourceforge.net/bwa.shtml))  at the same location as the reference! 

    [Short/Long]ReadRemap.sh] -b <BAM> -r <Modified_GRCh38_Ref> -d <bed_regions> -o <out_prefix>
    
 If CRAM file is used instead of BAM, then it is required to use the [original GRCh38 reference](https://bcm.box.com/s/ym4x3z61ib4okbguy7zn8lre8uo6mcxz) in the script with `-f` parameter.
 
    [Short/Long]ReadRemap.sh -b <CRAM> -r <Modified_GRCh38_Ref> -d <bed_regions> -f <Original_GRCh38_Ref> -o <out_prefix>
    
 To run with multiple threads, use `-t <no_threads>`
    
## Output

    <out_prefix>_original.sorted.bam <- Sorted BAM files corresponds to the extracted paired-end reads (original mapping)"
    <out_prefix>_remapped.sorted.bam <- Sorted bam file after BWA MEM step i.e. re-mapping of extracted reads to GRCh38 masked v2 reference (final output)"
