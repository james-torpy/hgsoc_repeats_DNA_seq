#!/bin/bash

# activate conda env using "conda activate general" before running!

home_dir="/share/ScratchGeneral/jamtor"
in_dir="$home_dir/projects/hgsoc_repeats/DNA-seq/results/bp_assoc_RT/tables"

# sort and index overall gtfs:
echo "Sorting and indexing gtfs in $in_dir..."

for f in $in_dir/*.gtf; do
  prefix=$(echo $f | sed "s/.gtf//")
  igvtools sort $f $prefix.sorted.gtf
  igvtools index $prefix.sorted.gtf
done;

# sort and index sample-specific gtfs:
for dir in $in_dir/AOCS*; do

  echo "Sorting and indexing gtfs in $dir..."

  for f in $dir/*.gtf; do

  	prefix=$(echo $f | sed "s/.gtf//")
  	igvtools sort $f $prefix.sorted.gtf

  	igvtools index $prefix.sorted.gtf

  done;

done;

rm $in_dir/*sorted.sorted*
rm $in_dir/**/*sorted.sorted*

