#!/bin/bash

##Written by Brian Bushnell
##Last modified April 8, 2020
##Description:  Summarizes all the runs of "processCorona.sh" in this directory.
##Then tars them if you want to send them somewhere.
##Should be run only after all samples are processed individually.

##Call variants in multisample mode to solve things.
##This reports the genotype of *all* samples at any position at which a variant is called in *any* sample.
callvariants.sh *_deduped_trimclip.sam.gz ref=NC_045512.fasta multisample out=allVars.vcf ow -Xmx4g usebias=f strandedcov minstrandratio=0 maf=0.6

##Make a summary of coverage at varous depth cutoffs for all libraries.
summarizecoverage.sh *basecov.txt out=coverageSummary.txt

mkdir output
cp *.sh output
cp *.bam* output
cp *.txt output
cp *.vcf output
cp *.fa output

rm results.tar
tar -cf results.tar output
