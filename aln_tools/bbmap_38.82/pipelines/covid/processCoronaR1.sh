#!/bin/bash

#Written by Brian Bushnell
#Last modified April 4, 2020
#Description:  Calls coronavirus variants from single-ended 300bp amplicon data.
#Usage:  processCoronaR1.sh <prefix>
#For example, "processCoronaR1.sh Sample1" if the data is in Sample1.fq.gz


NAME="$1"

#This line is in case the script is being re-run, clears the old junk
rm "$NAME"*.sam.gz "$NAME"*.bam "$NAME"*.bai "$NAME"*.txt "$NAME"*.fq.gz "$NAME"*.fa "$NAME"*.vcf  

#Recalibrate quality scores prior to any trimming (see recal.sh for details)
bbduk.sh in="$NAME".fastq.gz out="$NAME"_recal.fq.gz recalibrate -Xmx1g ow

#Perform adapter-trimming on the reads
#Also do quality trimming and filtering
bbduk.sh in="$NAME"_recal.fq.gz out="$NAME"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=truseq maq=14 qtrim=r trimq=10 maxns=0 ow -Xmx1g

#Align singleton reads to the reference
bbmap.sh ref=NC_045512.fasta in="$NAME"_trimmed.fq.gz outm="$NAME"_trimmed.sam.gz nodisk local maxindel=10 -Xmx4g

#Deduplicate based on mapping coordinates
dedupebymapping.sh in="$NAME"_trimmed.sam.gz out="$NAME"_deduped.sam.gz -Xmx4g

bbduk.sh in="$NAME"_deduped.sam.gz trimclip out="$NAME"_deduped_trimclip.sam.gz -Xmx1g

#Make bam/bai files for visualization
samtools view -bShu "$NAME"_deduped_trimclip.sam.gz | samtools sort -m 2G -@ 3 - -o "$NAME"_deduped_sorted.bam
samtools index "$NAME"_deduped_sorted.bam

#Call variants from the sam files
callvariants.sh in="$NAME"_deduped_trimclip.sam.gz ref=NC_045512.fasta out="$NAME"_vars.vcf -Xmx4g

#Generate a mutant reference by applying the detected variants to the reference
#This is essentially the reference-guided assembly of the strain
applyvariants.sh in=NC_045512.fasta out="$NAME"_genome.fa vcf="$NAME"_vars.vcf

#The remaining steps are just to allow visualization (for example, in IGV).

#At this point, "$NAME"_deduped_sorted.bam, "$NAME"_deduped_sorted.bam.bai, NC_045512.fasta, and "$NAME"_vars.vcf can be used for visualization in IGV.

#Optional step if you want to look at your own plot of coverage.
#Some regions of the virus, particularly the end, have low or zero coverage.
pileup.sh in="$NAME"_deduped_trimclip.sam.gz bincov="$NAME"_bincov.txt basecov="$NAME"_basecov.txt binsize=100 -Xmx4g

