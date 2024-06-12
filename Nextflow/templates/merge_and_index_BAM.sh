#!/bin/bash

# Script adapted from SNP simulation pipeline developed by L. Weber et al (2021) doi:10.1093/gigascience/giab062

# ------------------------------------------------
# Shell script to merge and index parsed BAM files
# ------------------------------------------------

# arguments:
# $1: input bams
# $2: merged bam
# $3  outdir

# -----------------------------------
# start runtime
# start=`date +%s`
# -----------------------------------


mkdir -p $3

# merge BAM files
samtools merge -@ 4 $3/$2 \
$1

# index merged BAM
samtools index $3/$2


# -----------------------------------
# end runtime
# end=`date +%s`
# runtime=`expr $end - $start`

# save runtime
# mkdir -p $1/merge_and_index_BAM
# echo runtime: $runtime seconds > $1/merge_and_index_BAM/runtime_merge_and_index_BAM.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
# mkdir -p $2/merge_and_index_BAM
# date > $2/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt
# -----------------------------------

