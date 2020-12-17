project_name="hgsoc_repeats"
exp_name="DNA-seq"

sample=$1

home_dir="/share/ScratchGeneral/jamtor/"
project_dir="$home_dir/projects/$project_name/$exp_name/"
manta_dir="results/manta/"
bwa_dir="results/bwa/"
genome_dir="genome/"

source /home/jamtor/.bashrc
conda activate py2.7

mkdir -p $manta_dir/AOCS-063-sub

configManta.py \
--region chr1 --region chr2 --region chr3 \
--region chr4 --region chr5 --region chr6 \
--region chr7 --region chr8 --region chr9 \
--region chr10 --region chr11 --region chr12 \
--region chr13 --region chr14 --region chr15 \
--region chr16 --region chr17 --region chr18 \
--region chr19 --region chr20 --region chr21 \
--region chr22 --region chrX --region chrY \
--normalBam ../../$bwa_dir/$sample\-1.bam \
--tumorBam ../../$bwa_dir/$sample\-5.bam \
--referenceFasta ../../$genome_dir/GRCh38.primary_assembly.genome.fa \
--runDir $project_dir/$manta_dir/$sample

../../$manta_dir/AOCS-063-sub/runWorkflow.py -m local -j 15