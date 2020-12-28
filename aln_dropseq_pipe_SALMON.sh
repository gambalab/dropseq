#!/bin/bash

# qsub -l nodes=1:ppn=15 -l naccesspolicy=singleuser -q fourjobs -F "sample_name path/to/read1 path/to/read2 number_of_target_cells project_out_folder" /home/users/ngs/dropseq_pipe/aln_dropseq_pipe_STARsolo.sh
#
# Notes:
# sample_name = name of the sample. It will be used to store results in $project_out_folder
# number_of_target_cells = usually equal to 2000 in Gaetano's experiments


####################
# Step 0: configs  #
####################

export PATH=/opt/software/java/jdk8/bin/:$PATH
wd="/home/tigem/gambardella/dropseq/" # path where the script is
cpus=15
xmx=96g

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

# genome and SALMON paths
export LD_LIBRARY_PATH=${wd}/aln_tools/salmon_14/lib:LD_LIBRARY_PATH
genomeAnn='/home/tigem/gambardella/Downloads/genome_idx/salmon_idx/salmon_index'
txp2gene='/home/tigem/gambardella/Downloads/genome_idx/salmon_idx/salmon_index/txp2gene.tsv'
salmon_aln="${wd}/aln_tools/salmon/bin/salmon"

# Paths to run emptyDrop
rscript_bin='/usr/bin/Rscript'
R_emptyDrop="${wd}/emptyDrop.R"


############################
# Step 0: paths and files  #
############################

echo "Analazyng sample = ${sample}"

mkdir -p "${outDir}/${sample}"
mkdir -p "${outDir}/${sample}/processed_fastq"
mkdir -p "${outDir}/${sample}/stats"
mkdir -p "${outDir}/${sample}/aln_alevin"
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

##########################################
# Step 2: Identify correct cell barcodes #
##########################################
groupBC  ${R1} ${sample} ${outDir}/${sample}/aln_alevin/BC_estimation
zcat ${outDir}/${sample}/aln_alevin/BC_estimation/${sample}_BC_counts.gz | awk -F " " '{ if($1>=5000) {print $2}}' > "${outDir}/${sample}/aln_alevin/BC_estimation/${sample}_whitelist.txt"

#####################
# Step 3: Map reads #
#####################

# Aligment and quantification with SALMON
#----------------------------------------
echo "Aln reads.."
ionice -c 3 ${salmon_aln} alevin -l ISR -1 ${R1} -2 ${R2} --dropseq -i ${genomeAnn} -p ${cpus} -o ${outDir}/${sample}/aln_alevin/ --tgMap ${txp2gene} --whitelist "${outDir}/${sample}/aln_alevin/BC_estimation/${sample}_whitelist.txt" --dumpMtx

# rm fastq
#ionice -c 3 rm ${R1} ${R2}

