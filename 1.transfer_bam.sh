#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/DNA-seq"
out_dir="$project_dir/raw_files/bams"

gadi_dir="/g/data1a/ku3/jt3341"
in_dir="jt3341/projects/hgsoc_repeats/DNA-seq/raw_files"
temp_dir="$gadi_dir/mdss/hgsoc_repeats/DNA-seq/raw_files"

sample_name=$1
tumour_only=$2

#sample_name="AOCS-063-sub"

mkdir -p $out_dir
cd $out_dir

# sleep for random time between 0 and 5 min to avoid Gadi rsync cap:
sl=$(( ( RANDOM % 300 )  + 1 ))
echo "Sleeping for $sl seconds to avoid rsync issue..."
sleep $sl

if [ $tumour_only = "true" ]; then

	filename=$sample_name\-1.bam.gz
  
    echo "Fetching $filename from mdss..."
    ssh jt3341@gadi.nci.org.au mdss get $in_dir/$filename $temp_dir
    
    echo "Transferring $filename from gadi..."
    rsync -avPS jt3341@gadi-dm.nci.org.au:$temp_dir/$filename $out_dir
    
    echo "Deleting $filename from gadi..."
    ssh jt3341@gadi.nci.org.au rm $temp_dir/$filename

else

  for i in 1 5; do 
  
    filename=$sample_name\-$i.bam.gz
  
    echo "Fetching $filename from mdss..."
    ssh jt3341@gadi.nci.org.au mdss get $in_dir/$filename $temp_dir
    
    echo "Transferring $filename from gadi..."
    rsync -avPS jt3341@gadi-dm.nci.org.au:$temp_dir/$filename $out_dir
    
    echo "Deleting $filename from gadi..."
    ssh jt3341@gadi.nci.org.au rm $temp_dir/$filename
    
  done

fi;

