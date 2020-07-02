#!/bin/bash

##Written by Brian Bushnell
##Last modified April 8, 2020
##Description:  Calls coronavirus variants from single-ended amplicon data.

##Usage:  processCoronaSE.sh <prefix>
##For example, "processCoronaSE.sh Sample1" if the data is in Sample1.fq.gz

##Set minimum coverage for genotype calls.
##Areas below this depth will be set to N in the consensus genome.
MINCOV=3

##Capture sample name from command line
NAME="$1"

##This line is in case the script is being re-run, to clears the old output
rm "$NAME"*.sam.gz "$NAME"*.bam "$NAME"*.bai "$NAME"*.txt "$NAME"*.fa "$NAME"*.vcf  

##Recalibrate quality scores prior to any trimming.
##Requires a recalibration matrix in the working directory (see recal.sh for details).
##This step is optional but useful for Illumina binned quality scores.
bbduk.sh in="$NAME".fastq.gz out="$NAME"_recal.fq.gz recalibrate -Xmx1g ow


##Perform adapter-trimming on the reads.
##Also do quality trimming and filtering.
##If desired, also do primer-trimming here by adding, e.g., 'ftl=20' to to trim the leftmost 20 bases.
bbduk.sh in="$NAME"_recal.fq.gz out="$NAME"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=truseq maq=14 qtrim=r trimq=10 maxns=0 ow -Xmx1g

##Align reads to the reference.
##NC_045512.fasta is the coronavirus reference, equivalent to bbmap/resources/Covid19_ref.fa
bbmap.sh ref=NC_045512.fasta in="$NAME"_trimmed.fq.gz outm="$NAME"_trimmed.sam.gz nodisk local maxindel=20 -Xmx4g ow

##Deduplicate based on mapping coordinates.
##Note that single-ended amplicon data will lose most of its data here.
dedupebymapping.sh in="$NAME"_trimmed.sam.gz out="$NAME"_deduped.sam.gz -Xmx4g ow

##Trim soft-clipped bases.
bbduk.sh in="$NAME"_deduped.sam.gz trimclip out="$NAME"_deduped_trimclip.sam.gz -Xmx1g ow

##Call variants from the sam files.
##The usebias=f/minstrandratio=0 flags *may be* necessary due to amplicon due to strand bias,
##and should be removed if the data is exclusively shotgun/metagenomic or otherwise randomly fragmented,
##and exhibits no strand bias.
callvariants.sh in="$NAME"_deduped_trimclip.sam.gz ref=NC_045512.fasta out="$NAME"_vars.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV"

##Calculate coverage.
pileup.sh in="$NAME"_deduped_trimclip.sam.gz bincov="$NAME"_bincov.txt basecov="$NAME"_basecov.txt binsize=100 -Xmx4g ow

##Calculate reduced coverage as per CallVariants defaults (ignoring outermost 5bp of reads).
pileup.sh in="$NAME"_deduped_trimclip.sam.gz basecov="$NAME"_basecov_border5.txt -Xmx4g ow border=5

##Generate a mutant reference by applying the detected variants to the reference.
##This is essentially the reference-guided assembly of the strain.
##Also changes anything below depth MINCOV to N (via the mindepth flag).
applyvariants.sh in=NC_045512.fasta out="$NAME"_genome.fa vcf="$NAME"_vars.vcf basecov="$NAME"_basecov.txt ow mindepth="$MINCOV"

##Make bam/bai files; requires samtools to be installed.
##This step is only necessary for visualization, not variant-calling.
samtools view -bShu "$NAME"_deduped_trimclip.sam.gz | samtools sort -m 2G -@ 3 - -o "$NAME"_deduped_sorted.bam
samtools index "$NAME"_deduped_sorted.bam

##At this point, "$NAME"_deduped_sorted.bam, "$NAME"_deduped_sorted.bam.bai, NC_045512.fasta, and "$NAME"_vars.vcf can be used for visualization in IGV.
