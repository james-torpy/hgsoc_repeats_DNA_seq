#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
genome_dir="$home_dir/projects/hgsoc_repeats/DNA-seq/genome"

genes=( "MAML2" "EXT2" "FGF18" "TOX" "UBE2C" "POLD1" "EPHA3" "FCRL4" "ANKRD11" "ATM" "EML4" "MECOM" "CXCR2" )
annot_file="gencode.v35.basic.annotation.gtf"
out_file="custom.genes.gencode.v35.gtf"

for g in ${genes[@]}; do
  grep $g $genome_dir/$annot_file >> $genome_dir/$out_file
done;
