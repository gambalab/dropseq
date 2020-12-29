#!/bin/bash

# qsub -l nodes=1:ppn=15 -l naccesspolicy=singleuser -q fourjobs -F "sample_name path/to/read1 path/to/read2 number_of_target_cells project_out_folder" /home/users/ngs/dropseq_pipe/aln_dropseq_pipe_STARsolo.sh
#
# Notes:
# sample_name = name of the sample. It will be used to store results in $project_out_folder
# number_of_target_cells = usually equal to 2000 in Gaetano's experiments


####################
# Step 0: configs  #
####################

#export PATH=/opt/software/java/jdk8/bin/:$PATH
wd="/home/tigem/gambardella/dropseq_pipe/" # path where the script is
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

# UMI tool
umitool="/usr/local/bin/umi_tools"

# genome and STAR paths
genomeAnn='/home/tigem/gambardella/Downloads/genome_idx/human_cellecta_myco'
gtfAnn='/home/tigem/gambardella/Downloads/genome_idx/gencode.v36.primary_assembly.annotation.gtf'
star="${wd}/aln_tools/STAR-2.7.5a/STAR"

# subread
fcount="${wd}/aln_tools/subread/bin/featureCounts"

############################
# Step 0: paths and files  #
############################

echo "Analazyng sample = ${sample}"

mkdir -p "${outDir}/${sample}"
mkdir -p "${outDir}/${sample}/processed_fastq"
mkdir -p "${outDir}/${sample}/stats"
mkdir -p "${outDir}/${sample}/aln_fcount"
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
groupBC  ${R1} ${sample} ${outDir}/${sample}/aln_fcount/BC_estimation
zcat ${outDir}/${sample}/aln_fcount/BC_estimation/${sample}_BC_counts.gz | awk -F " " '{ if($1>=1000) {print $2}}' > "${outDir}/${sample}/aln_fcount/BC_estimation/${sample}_whitelist.txt"

###########################################################
# Step 3: Extract barcdoes and UMIs and add to read names #
###########################################################

ionice -c 3 ${umitool} extract --bc-pattern=CCCCCCCCCCCCNNNNNNNN \
                  --stdin ${R1} \
                  --stdout "${outDir}/${sample}/processed_fastq/${sample}_R1_extracted_BC.fastq.gz" \
                  --read2-in ${R2} \
                  --read2-out="${outDir}/${sample}/processed_fastq/${sample}_R2_extracted_BC.fastq.gz" \
                  --filter-cell-barcode \
                  --whitelist="${outDir}/${sample}/aln_fcount/BC_estimation/${sample}_whitelist.txt";

# rm fastq
#ionice -c 3 rm ${R1} ${R2}


#####################
# Step 4: Map reads #
#####################
ionice -c 3 ${star} --runThreadN ${cpus} \
     --genomeDir "${genomeAnn}" \
     --readFilesIn "${outDir}/${sample}/processed_fastq/${sample}_R2_extracted_BC.fastq.gz" \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outFileNamePrefix "${outDir}/${sample}/aln_fcount/" \
     --outSAMtype BAM SortedByCoordinate;

#################################
# Step 5: Assign reads to genes #
#################################
ionice -c 3 ${fcount} -a "${gtfAnn}" \
              -o "${outDir}/${sample}/aln_fcount/gene_assigned" \
              -R BAM "${outDir}/${sample}/aln_fcount/Aligned.sortedByCoord.out.bam" \
              -T ${cpus};


# sort BAM
echo "Sort BAM file and Index sorted BAM"
ionice -c 3 rm "${outDir}/${sample}/aln_fcount/Aligned.sortedByCoord.out.bam"
ionice -c 3 samtools sort --threads ${cpus} "${outDir}/${sample}/aln_fcount/Aligned.sortedByCoord.out.bam.featureCounts.bam" -o "${outDir}/${sample}/aln_fcount/${sample}.featureCounts.srt.bam"
ionice -c 3 samtools index -@ ${cpus} "${outDir}/${sample}/aln_fcount/${sample}.featureCounts.srt.bam"
ionice -c 3 rm "${outDir}/${sample}/aln_fcount/Aligned.sortedByCoord.out.bam.featureCounts.bam"

########################################
# Step 6: Count UMIs per gene per cell #
########################################
ionice -c 3 ${umitool} count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I "${outDir}/${sample}/aln_fcount/${sample}.featureCounts.srt.bam" -S "${outDir}/${sample}/aln_fcount/${sample}.counts.tsv.gz" -L "${outDir}/${sample}/aln_fcount/${sample}.umi.count.txt"



