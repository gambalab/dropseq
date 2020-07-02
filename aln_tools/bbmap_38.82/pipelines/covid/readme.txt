Written by Brian Bushnell
Last modified April 8, 2020

Contents:
This directory contains a collection of scripts written for calling variants and generating consensus genomes from SARS-CoV-2 ("Covid") Illumina data using BBTools.  They were designed and optimized for several sets of libraries (from NextSeq and NovaSeq platforms) containing a mix of shotgun and amplicon data using various primers, and may need adjustment for your specific experimental design.
These scripts, and those in BBTools, should work fine in Bash on Linux, MacOS, or Windows 10.  Other shells generally work fine.

Usage:
These scripts are just guidelines showing how I am processing Covid data.
If you want to use them without modification, follow these steps:

1) Install BBTools (you need the latest version - 38.81+ - due to some new features I added for Covid!).
   a) If you don't have Java, install Java.
      i) If you run java --version on the command line and it reports version 8+, you're fine.
      ii) Otherwise, you can get it from https://openjdk.java.net/install/index.html
   b) Download BBTools 38.81 or higher from https://sourceforge.net/projects/bbmap/
   c) Unzip the archive: tar -xzf BBMap_38.81.tar.gz
   d) Add it to the path: export PATH=/path/to/bbmap/:$PATH
      i) Now you can run the shell scripts from wherever.
   e) Samtools is not necessary, but recommended if you want to make the bam files for visualization.
2) Rename and interleave the files if they are paired, so libraries are in the format "prefix.fq.gz" with one file per library.
3) Copy the Covid reference to the directory.
4) Modify the template processCoronaWrapper.sh:
   a) Pick a library to use for quality-score calibration if desired (line 17).
   b) Add lines, or make a loop, so that processCorona.sh is called on all of the libraries (lines 21-22).
   c) Delete lines 8-9 to let the script run.
5) Run processCoronaWrapper.sh.

Note:
You can remove human reads with a command like this (where reads.fq can be single-ended or paired/interleaved):
bbmap.sh ref=hg19.fa in=reads.fq outm=human.fq outu=nonhuman.fq bloom

I don't provide a script for that because the data I am processing has already had human reads removed.

