#!/bin/bash

# need to conda activate general:
source ~/.bashrc
conda activate general

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/DNA-seq"
in_dir="$project_dir/results/bwa"
genome_dir="$project_dir/genome"

#######
#script_dir="$project_dir/scripts"
#log_dir="$project_dir/logs/calculate_bam_coverage"
#mkdir -p $log_dir
#
#for f in $in_dir/*.bam; do
#
#	sample_name=$(echo $f | sed "s/.bam//" | sed "s/^.*\///")
#
#	if [ ! -f "$in_dir/$sample_name.coverage.bed" ]; then
#
#	  echo "Calculating coverage of $f.."
#
#	  qsub -pe smp 6 -N $sample_name.genomecov -wd $log_dir -b y -j y -V \
#	    -P TumourProgression $script_dir/calculate_coverage_from_bam.sh $sample_name
#
#	else
#
#	  echo "Coverage already calculated for $f..."
#
#	fi;
#
#done;
######

sample_name=$1

echo "Calculating coverage from $sample_name.bam..."

# calculate coverage across bam file:
time bedtools genomecov -bga -ibam $in_dir/$sample_name.bam > \
  $in_dir/$sample_name.coverage.bed
 echo "Calculating overall mean coverage from $prefix.coverage.bed..."

# calculate overall mean coverage:
time awk '{total+=$4} END {print total/NR}' $in_dir/$sample_name.coverage.bed > \
  $in_dir/$sample_name.mean.coverage.txt

#  # sort bed file for igv:
#  sort -k 1,1 -k 2,2n -k 3,3n $in_dir/$prefix.coverage.bed | bgzip -c > \
#    $in_dir/$prefix.coverage.bed.gz
#  
#  # and index:
#  tabix -pbed $in_dir/$prefix.coverage.bed.gz 


