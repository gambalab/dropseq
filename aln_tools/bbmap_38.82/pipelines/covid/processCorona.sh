#!/bin/bash

##Written by Brian Bushnell
##Last modified April 17, 2020
##Description:  Calls coronavirus variants from interleaved PE amplicon data.
##This script assumes input data is paired.

##Usage:  processCorona.sh <prefix>
##For example, "processCorona.sh Sample1" if the data is in Sample1.fq.gz

##Set minimum coverage for genotype calls.
##Areas below this depth will be set to N in the consensus genome.
MINCOV=3

##Grab the sample name from the command line
NAME="$1"

##This line is in case the script is being re-run, to clear the old output
rm "$NAME"*.sam.gz "$NAME"*.bam "$NAME"*.bai "$NAME"*.txt "$NAME"*.fa "$NAME"*.vcf

##If data is paired in twin files, interleave it into a single file.
##Otherwise, skip this step.
##In this case, the files are assumed to be named "Sample1_R1.fq.gz" and "Sample1_R2.fq.gz"
#reformat.sh in="$NAME"_R#.fq.gz out="$NAME".fq.gz

##Recalibrate quality scores prior to any trimming.
##Requires a recalibration matrix in the working directory (see recal.sh for details).
##This step is optional but useful for Illumina binned quality scores.
bbduk.sh in="$NAME".fq.gz out="$NAME"_recal.fq.gz recalibrate -Xmx1g ow

##Discover adapter sequence for this library based on read overlap.
##You can examine the adapters output file afterward if desired;
##If there were too few short-insert pairs this step will fail (and you can just use the default Illumina adapters).
bbmerge.sh in="$NAME"_recal.fq.gz outa="$NAME"_adapters.fa ow reads=1m

##Remove duplicates by sequence similarity.
##This is more memory-efficient than dedupebymapping.
clumpify.sh in="$NAME"_recal.fq.gz out="$NAME"_clumped.fq.gz zl=9 dedupe s=2 passes=4 -Xmx31g

##Perform adapter-trimming on the reads.
##Also do quality trimming and filtering.
##If desired, also do primer-trimming here by adding, e.g., 'ftl=20' to to trim the leftmost 20 bases.
##If the prior adapter-detection step failed, use "ref=adapters"
bbduk.sh in="$NAME"_clumped.fq.gz out="$NAME"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref="$NAME"_adapters.fa maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx1g

##Split into Covid and non-Covid reads if this has not already been done.
##This step can be skipped if non-Covid was already removed.
##NC_045512.fasta is the coronavirus reference, equivalent to bbmap/resources/Covid19_ref.fa
bbduk.sh ow -Xmx1g in="$NAME"_trimmed.fq.gz ref=NC_045512.fasta outm="$NAME"_covid.fq.gz outu="$NAME"_noncovid.fq.gz k=25

##Align reads to the reference.
bbmap.sh ref=NC_045512.fasta in="$NAME"_trimmed.fq.gz outm="$NAME"_mapped.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12

##Deduplicate based on mapping coordinates.
##Note that if you use single-ended amplicon data, you will lose most of your data here.
dedupebymapping.sh in="$NAME"_mapped.sam.gz out="$NAME"_deduped.sam.gz -Xmx31g ow

#Remove junk reads with unsupported unique deletions; these are often chimeric.
filtersam.sh ref=NC_045512.fasta ow in="$NAME"_deduped.sam.gz out="$NAME"_filtered.sam.gz mbad=1 del sub=f mbv=0 -Xmx4g

#Remove junk reads with multiple unsupported unique substitutions; these are often junk, particularly on Novaseq.
#This step is not essential but reduces noise.
filtersam.sh ref=NC_045512.fasta ow in="$NAME"_filtered.sam.gz out="$NAME"_filtered2.sam.gz mbad=1 sub mbv=2 -Xmx4g

##Trim soft-clipped bases.
bbduk.sh in="$NAME"_filtered2.sam.gz trimclip out="$NAME"_trimclip.sam.gz -Xmx1g ow

##Call variants from the sam files.
##The usebias=f/minstrandratio=0 flags are necessary due to amplicon due to strand bias,
##and should be removed if the data is exclusively shotgun/metagenomic or otherwise randomly fragmented.
callvariants.sh in="$NAME"_trimclip.sam.gz ref=NC_045512.fasta out="$NAME"_vars.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV" minedistmax=30 minedist=16

##Calculate coverage.
pileup.sh in="$NAME"_trimclip.sam.gz bincov="$NAME"_bincov.txt basecov="$NAME"_basecov.txt binsize=100 -Xmx4g ow

##Calculate reduced coverage as per CallVariants defaults (ignoring outermost 5bp of reads).
pileup.sh in="$NAME"_trimclip.sam.gz basecov="$NAME"_basecov_border5.txt -Xmx4g ow border=5

##Generate a mutant reference by applying the detected variants to the reference.
##This is essentially the reference-guided assembly of the strain.
##Also changes anything below depth MINCOV to N (via the mindepth flag).
applyvariants.sh in=NC_045512.fasta out="$NAME"_genome.fa vcf="$NAME"_vars.vcf basecov="$NAME"_basecov.txt ow mindepth="$MINCOV"

##Make bam/bai files; requires samtools to be installed.
##This step is only necessary for visualization, not variant-calling.
samtools view -bShu "$NAME"_trimclip.sam.gz | samtools sort -m 2G -@ 3 - -o "$NAME"_sorted.bam
samtools index "$NAME"_sorted.bam

##At this point, "$NAME"_sorted.bam, "$NAME"_sorted.bam.bai, NC_045512.fasta, and "$NAME"_vars.vcf can be used for visualization in IGV.

