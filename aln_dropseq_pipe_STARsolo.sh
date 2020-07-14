#!/bin/sh
## Set the working Directorty to be the current one, i.e. from where you are submitting the script
#$ -cwd

# qsub -l nodes=1:ppn=15 -l naccesspolicy=singleuser -q fourjobs -F "sample_name path/to/read1 path/to/read2 number_of_target_cells project_out_folder" /home/users/ngs/dropseq_pipe/drop_seq/aln_dropseq_pipe_STARsolo.sh
#
# Notes:
# sample_name = name of the sample. It will be used to store results in $project_out_folder
# number_of_target_cells = usually equal to 2000 in Gaetano's experiments


####################
# Step 0: configs  #
####################

export PATH=/opt/software/java/jdk8/bin/:$PATH
wd="/home/users/ngs/dropseq_pipe" # path where the script is
cpus=15
xmx=48g

sample=$1
read1=$2
read2=$3
ncells=$4
outDir=$5


####################
# Tools            #
####################

# Dropseq lib path
. ${wd}/dropseq_aln_v2.lib

# dropseq tool
dropseq_tool="${wd}/aln_tools/Drop-seq_tools-2.3.0/jar/dropseq.jar"

# Picard tools
picard_tool="${wd}/aln_tools/picard/picard.jar"

# bbduk
bbmap="${wd}/aln_tools/bbmap_38.82/" # <-- Put / at the end
bbmapReource="${wd}/aln_tools/bbmap_38.82/resources"

# genome and STAR paths
genomeAnn='/home/users/ngs/references/dropseq/dropseq_genome_idx'
gtfAnn='/home/users/ngs/references/dropseq/genome_annotations/gencode.v34.primary_assembly.annotation.gtf'
star='/opt/software/ngs/bin/STAR275a'

# Paths to run emptyDrop
rscript_bin='/opt/software/R/stable/3.5.1/bin/Rscript'
R_emptyDrop='/home/users/gambardella/work/exome_variants/annovarAnnotate.R'


############################
# Step 0: paths and files  #
############################

echo "Analazyng sample = ${sample}"

mkdir -p "${outDir}/${sample}"
mkdir -p "${outDir}/${sample}/processed_fastq"
mkdir -p "${outDir}/${sample}/stats"
mkdir -p "${outDir}/${sample}/aln"
cd "${outDir}/${sample}"


############################
# Step 1: Reads filtering  #
############################

# FASTQ to Bam 
# -------------
fastq2sam ${xmx} ${read1} ${read2} ${sample} ${outDir}/${sample}/processed_fastq


# TAG low quality BARCODES 
# -------------------------
# mark barcodes with at least one base with a quality less then 10
tagLowQualityReads 'bc' '1-12' '10' ${xmx} ${outDir}/${sample}/processed_fastq ${outDir}/${sample}/stats ${sample}
ionice -c 3 rm "${outDir}/${sample}/processed_fastq/${sample}_unl_read_pairs.bam"

# TAG low quality UMI 
#---------------------
# mark UMI with at least one base with a quality less then 10
tagLowQualityReads 'umi' '13-20' '10' ${xmx} ${outDir}/${sample}/processed_fastq/ ${outDir}/${sample}/stats/ ${sample}
ionice -c 3 rm "${outDir}/${sample}/processed_fastq/${sample}_unl_read_pairs_tagged_cells.bam"


# Filter Out Low quality barcodes and UMI 
# ----------------------------------------
# discard low quality marked UMI and Barcodes
filterBAM ${xmx} ${outDir}/${sample}/processed_fastq/ ${sample}
ionice -c 3 rm "${outDir}/${sample}/processed_fastq/${sample}_unl_read_pairs_tagged_cells_umi.bam"


# Trim adapter from 5' of reads 2 
#---------------------------------
trimReads 'adapter' ${xmx} ${outDir}/${sample}/processed_fastq/ ${outDir}/${sample}/stats/ ${sample}
ionice -c 3 rm "${outDir}/${sample}/processed_fastq/${sample}_unl_read_pairs_filtered.bam"


# Trim polyA from 3' of reads 2 
#-------------------------------
trimReads 'polyA' ${xmx} ${outDir}/${sample}/processed_fastq/ ${outDir}/${sample}/stats/ ${sample}
ionice -c 3 rm "${outDir}/${sample}/processed_fastq/${sample}_unl_read_pairs_filtered_adapter_trimmed.bam"


# convert everythimg to FASTQ again 
#-----------------------------------
sam2fastq ${xmx} ${outDir}/${sample}/processed_fastq/ ${sample}
ionice -c 3 rm "${outDir}/${sample}/processed_fastq/${sample}_unl_read_pairs_filtered_adapter_polyA_trimmed.bam"

# Get survived Barcodes from R2 
# ------------------------------
intersectFASTQ ${read1} 20000000 ${xmx} ${outDir}/${sample}/processed_fastq/ ${sample}

R1="${outDir}/${sample}/processed_fastq/${sample}_filtered_and_trimmed_R1.fastq.gz"
R2="${outDir}/${sample}/processed_fastq/${sample}_filtered_and_trimmed_R2.fastq.gz"

# Aligment and quantification with STARsolo 
# ------------------------------------------
echo "Aln reads.."
ionice -c 3 ${star} --soloType CB_UMI_Simple \
     --soloCBstart 1 \
     --soloCBlen 12 \
     --soloUMIstart 13 \
     --soloUMIlen 8 \
     --soloUMIdedup 1MM_Directional \
     --soloCBwhitelist None \
     --soloFeatures Gene Velocyto \
     --runThreadN ${cpus} \
     --genomeDir ${genomeAnn} \
     --readFilesIn ${R2} ${R1} \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outFileNamePrefix ${outDir}/${sample}/aln/ \
     --soloCellFilter CellRanger2.2 ${ncells} 0.99 10 \
     --outSAMtype BAM Unsorted;

ionice -c 3 ${R1} ${R2}


# Run emptyDrop with R 
# ------------------------------
echo "Filter out cell matrix with emptyDrop.."
${rscript_bin} --vanilla ${R_emptyDrop} "${outDir}/${sample}/aln/Solo.out/" ${wd} ${sample} ${cpus}

# Zip files 
# ------------------------------
echo "Zip files.."
gzip ${outDir}/${sample}/aln/Solo.out/Gene/raw/matrix.mtx
gzip ${outDir}/${sample}/aln/Solo.out/Gene/raw/barcodes.tsv
gzip ${outDir}/${sample}/aln/Solo.out/Gene/raw/features.tsv

gzip ${outDir}/${sample}/aln/Solo.out/Gene/filtered/matrix.mtx
gzip ${outDir}/${sample}/aln/Solo.out/Gene/filtered/barcodes.tsv
gzip ${outDir}/${sample}/aln/Solo.out/Gene/filtered/features.tsv
