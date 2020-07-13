#!/bin/sh

####################
# Step 0: configs  #
####################

#export PATH=/opt/software/java/jdk8/bin/:$PATH

sample='TS1'
cpus=6
xmx=48g
read1="/mnt/SSD_WIN10/fastq/TS1_R1_001.fastq.gz"
read2="/mnt/SSD_WIN10/fastq/TS1_R2_001.fastq.gz"
outDir="/mnt/SSD500/"
ncells=2000

####################
# Tools            #
####################

# Dropseq lib path
. /home/gg/work/package/drop_seq/dropseq_aln_v2.lib

# dropseq tool
dropseq_tool="./Drop-seq_tools-2.3.0/jar/dropseq.jar"

# Picard tools
picard_tool="./aln_tools/picard/picard.jar"

# bbduk
bbmap="./aln_tools/bbmap_38.82/" # <-- Put / at the end
bbmapReource="./aln_tools/bbmap_38.82/resources"

# genome and STAR paths
genomeAnn='/path/to/genome_STAR_idx_folder'
gtfAnn='/path/to/gtf'
star='./STAR-2.7.5a/STAR'

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
