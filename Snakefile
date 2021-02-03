# Run command:
# snakemake --reason --use-conda --cores 210 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N cnv_ident.smk -wd '/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/DNA-seq/logs' -b y -j y -V -P TumourProgression' -j 6

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script remaps bams to hg38 and uses SvABA and Manta to identify 
# breakpoints in HGSOC genomes ###

# create quota record file:

# check max quota used in quota_track.txt:
# awk '{print $1}' /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/DNA-seq/quota_track.txt | \
# grep -v GB | grep -v [a-zA-Z] | sort | uniq | tail -1

# max genomes than can be run at one one time (based on 10 TB quota): 6, and 5 in queue

# define directories:
project_name = 'hgsoc_repeats'
exp_name = 'DNA-seq'

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/' + exp_name + '/'
results_dir = project_dir + 'results/'
genome_dir = project_dir + 'genome/'
script_dir = project_dir + 'scripts/'


conda_dir = "/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/"
env_dir = conda_dir + "envs/py37/envs/snk/bin/"

bam_dir = 'raw_files/bams/'
temp_dir = project_dir + bam_dir + '/temp/'
fq_dir = 'raw_files/fastq/'
bwa_dir = 'results/bwa/'
svaba_dir = 'results/svaba/'
manta_dir = 'results/manta/'
manta_bin = '/g/data1a/ku3/jt3341/local/lib/manta-1.5.0/bin/'

#SAMPLES = list([
##    "AOCS-063-sub"
#	"AOCS-083", "AOCS-085", "AOCS-090", "AOCS-092",
#    "AOCS-063", "AOCS-064", "AOCS-065", "AOCS-075", 
#    "AOCS-076", "AOCS-077", "AOCS-078", "AOCS-080", 
#    "AOCS-083", "AOCS-084", "AOCS-085","AOCS-086", 
#    "AOCS-090", "AOCS-091", "AOCS-092", "AOCS-093", 
#    "AOCS-094", "AOCS-095", "AOCS-107", "AOCS-108", 
#    "AOCS-109", "AOCS-111", "AOCS-112", "AOCS-113", 
#    "AOCS-114", "AOCS-115", "AOCS-116", "AOCS-122", 
#    "AOCS-123", "AOCS-124", "AOCS-125", "AOCS-126", 
#    "AOCS-128", "AOCS-130", "AOCS-131", "AOCS-133", 
#    "AOCS-137"
#    #, "AOCS-143", "AOCS-144", "AOCS-145", 
##    "AOCS-146", "AOCS-147", "AOCS-148", "AOCS-149", 
##    "AOCS-153"
#])

# regenerate bams for:
SAMPLES = list([
  "AOCS-064", "AOCS-075", "AOCS-076", "AOCS-080", 
  "AOCS-083", "AOCS-085"
#  ,"AOCS-090", "AOCS-107", 
#  "AOCS-112", "AOCS-114", "AOCS-116",  "AOCS-122", 
])

## failed to transfer:
#SAMPLES = list([
#    "AOCS-152"
#])

#rule all:
#    input:
#        expand(
#            'completed_files/{sample}_bams_removed',
#            sample=SAMPLES
#        )

rule all:
    input:
        expand(
            bwa_dir + '{sample}-5.bam',
            bwa_dir + '{sample}-5.bam.bai',
            bwa_dir + '{sample}-1.bam',
            bwa_dir + '{sample}-1.bam.bai',
            sample=SAMPLES
        )

rule transfer:
    output:
        bam1 = 'raw_files/bams/{sample}-5.bam.gz',
        bam2 = 'raw_files/bams/{sample}-1.bam.gz'
    threads: 2
    shell:
        'mkdir -p logs/transfer; ' + 
        'cd logs/transfer; ' + 
        ' ../../scripts/1.transfer_bam.sh' + 
        ' {wildcards.sample}' +
            ' 2> {wildcards.sample}.transfer.errors'

rule gunzip:
    input:
        bamin1 = bam_dir + '{sample}-5.bam.gz',
        bamin2 = bam_dir + '{sample}-1.bam.gz'
    output:
        bamout1 = bam_dir + '{sample}-5.bam',
        bamout2 = bam_dir + '{sample}-1.bam'
    threads: 8
    shell:
        'mkdir -p logs/gunzip; ' + 
        'cd logs/gunzip; ' + 
        'pigz -d ../../{input.bamin1}' +
            ' 2> {wildcards.sample}.gunzip1.errors; ' + 
        env_dir + 'samtools quickcheck -v {output.bamout1} > ' + bam_dir + 
            '/{wildcards.sample}_bad_bams.fofn && echo "bam ok" || ' + 
            'echo "some files failed check, see bad_bams.fofn"; ' +
        'pigz -d ../../{input.bamin2}' +
            ' 2> {wildcards.sample}.gunzip1.errors; ' + 
        env_dir + 'samtools quickcheck -v {output.bamout2} > ' + bam_dir + 
            '/{wildcards.sample}_bad_bams.fofn && echo "bam ok" || ' + 
            'echo "some files failed check, see bad_bams.fofn"'

rule sort1:
    input:
        bamin1 = bam_dir + '{sample}-5.bam',
        bamin2 = bam_dir + '{sample}-1.bam',
    output:
        bamout1 = bam_dir + '{sample}-5.sorted.by.name.bam',
        bamout2 = bam_dir + '{sample}-1.sorted.by.name.bam'
    threads: 10
    shell:
        'mkdir -p logs/sort1; ' + 
        'cd logs/sort1; ' + 
        env_dir + 'samtools sort -n -@ 10 ../../{input.bamin1} -o ' + 
        	'../../{output.bamout1} 2> {wildcards.sample}.sort1a.errors; ' +
        'rm ../../{input.bamin1}; ' +
        env_dir + 'samtools sort -n -@ 10 ../../{input.bamin2} -o ' + 
        	'../../{output.bamout2} 2> {wildcards.sample}.sort1a.errors; ' +
        'rm ../../{input.bamin2}'

rule bam2fastq:
    input:
        bamin1 = bam_dir + '{sample}-5.sorted.by.name.bam',
        bamin2 = bam_dir + '{sample}-1.sorted.by.name.bam'
    output:
        fq1a = fq_dir + '{sample}-5-R1.fq',
        fq1b = fq_dir + '{sample}-5-R2.fq',
        fq_1single = fq_dir + '{sample}-5-singles.fq',
        fq2a = fq_dir + '{sample}-1-R1.fq',
        fq2b = fq_dir + '{sample}-1-R2.fq',
        fq_2single = fq_dir + '{sample}-1-singles.fq',
    threads: 15
    shell:
        'mkdir -p logs/bam2fastq; ' + 
        'cd logs/bam2fastq; ' + 
        env_dir + 'samtools bam2fq -1 ../../{output.fq1a} ' +
        ' -2 ../../{output.fq1b} -@ 15 -s ../../{output.fq_1single} ' +
        '../../{input.bamin1} 2> {wildcards.sample}.bam2fastq1.errors; ' + 
        'rm ../../{input.bamin1}; ' + 
        env_dir + 'samtools bam2fq -1 ../../{output.fq2a} ' +
        ' -2 ../../{output.fq2b} -@ 15 -s ../../{output.fq_2single} ' +
        '../../{input.bamin2} 2> {wildcards.sample}.bam2fastq1.errors; ' + 
        'rm ../../{input.bamin2}'

rule bwa_align:
    input:
        fq1a = fq_dir + '{sample}-5-R1.fq',
        fq1b = fq_dir + '{sample}-5-R2.fq',
        fq_1single = fq_dir + '{sample}-5-singles.fq',
        fq2a = fq_dir + '{sample}-1-R1.fq',
        fq2b = fq_dir + '{sample}-1-R2.fq',
        fq_2single = fq_dir + '{sample}-1-singles.fq',
        ind = genome_dir + 'GRCh38.primary_assembly.genome.fa.amb'
    output:
        sam1 = bwa_dir + '{sample}-5.sam',
        sam2 = bwa_dir + '{sample}-1.sam'
    threads: 35
    shell:
        'mkdir -p logs/bwa_align; ' + 
        'cd logs/bwa_align; ' + 
        env_dir + 'bwa mem -t 35 ' + genome_dir + 
            'GRCh38.primary_assembly.genome.fa ../../{input.fq1a} ' +
            '../../{input.fq1b} > ../../{output.sam1}; ' +
       	'rm ../../{input.fq1a}; ' + 
       	'rm ../../{input.fq1b}; ' + 
       	'rm ../../{input.fq_1single}; ' + 
        env_dir + 'bwa mem -t 35 ' + genome_dir + 
            'GRCh38.primary_assembly.genome.fa ../../{input.fq2a} ' +
            '../../{input.fq2b} > ../../{output.sam2}; ' +
       	'rm ../../{input.fq2a}; ' + 
       	'rm ../../{input.fq2b}; ' + 
       	'rm ../../{input.fq_2single}'

rule bam1:
    input:
        sam1 = bwa_dir + '{sample}-5.sam',
        sam2 = bwa_dir + '{sample}-1.sam'
    output:
        bam1 = bwa_dir + '{sample}-5.unsorted.bam',
        bam2 = bwa_dir + '{sample}-1.unsorted.bam'
    threads: 15
    shell:
        'mkdir -p logs/bam; ' + 
        'cd logs/bam; ' + 
        env_dir + 'samtools view -bh  -@ 15 ../../{input.sam1} > ' +
        	'../../{output.bam1}; ' +
        'rm ../../{input.sam1}; ' +
        env_dir + 'samtools view -bh  -@ 15 ../../{input.sam2} > ' +
        	'../../{output.bam2}; ' +
        'rm ../../{input.sam2}; '

rule sort2:
    input:
        unbam1 = bwa_dir + '{sample}-5.unsorted.bam',
        unbam2 = bwa_dir + '{sample}-1.unsorted.bam'
    output:
        bam1 = bwa_dir + '{sample}-5.bam',
        bam2 = bwa_dir + '{sample}-1.bam'
    threads: 10
    shell:
        'mkdir -p logs/sort2; ' + 
        'cd logs/sort2; ' + 
        env_dir + 'samtools sort -@ 10 ../../{input.unbam1} ' + 
        	'-o ../../{output.bam1} 2> {wildcards.sample}.sort2.errors; ' +
        'rm ../../{input.unbam1}; ' +
        env_dir + 'samtools sort -@ 10 ../../{input.unbam2} ' + 
        	'-o ../../{output.bam2} 2> {wildcards.sample}.sort2.errors; ' +
        'rm ../../{input.unbam2}'

rule index:
    input:
        bam1 = bwa_dir + '{sample}-5.bam',
        bam2 = bwa_dir + '{sample}-1.bam'
    output:
        bai1 = bwa_dir + '{sample}-5.bam.bai',
        bai2 = bwa_dir + '{sample}-1.bam.bai'
    threads: 15
    shell:
        'mkdir -p logs/index; ' + 
        'cd logs/index; ' + 
        env_dir + 'samtools index -@ 15 ../../{input.bam1}' +
            ' 2> {wildcards.sample}.index1.errors; ' + 
        env_dir + 'samtools index -@ 15 ../../{input.bam2}' +
            ' 2> {wildcards.sample}.index1.errors'

rule svaba:
   input:
       gbam = bwa_dir + '{sample}-5.bam',
       gbai = bwa_dir + '{sample}-5.bam.bai',
       sbam = bwa_dir + '{sample}-1.bam',
       sbai = bwa_dir + '{sample}-1.bam.bai'
   output:
       svaba_dir + '{sample}/{sample}.discordant.txt.gz'
   threads: 15
   shell:
       'mkdir -p logs/svaba; ' + 
        'cd logs/svaba; ' + 
        'mkdir -p ../../' + svaba_dir + '{wildcards.sample}/; ' +
        'svaba run -t ' + project_dir + '{input.sbam} -n ' + 
        project_dir + '{input.gbam} -G ' + 
            genome_dir + 'GRCh38.primary_assembly.genome.fa -a ../../' + 
            svaba_dir + 
            '{wildcards.sample}/{wildcards.sample} ' + 
            '-p 15' +
            ' 2> {wildcards.sample}.svaba.errors'

rule manta:
   input:
       gbam = bwa_dir + '{sample}-5.bam',
       gbai = bwa_dir + '{sample}-5.bam.bai',
       sbam = bwa_dir + '{sample}-1.bam',
       sbai = bwa_dir + '{sample}-1.bam.bai'
   output:
       manta_dir + '{sample}/results/variants/somaticSV.vcf.gz'
   threads: 15
   shell:
       'mkdir -p logs/manta; ' + 
        'cd logs/manta; ' + 
        '../../scripts/2.manta.sh ' + 
        '{wildcards.sample} ' +
            ' 2> {wildcards.sample}.manta.errors'

rule donefile:
   input:
        svaba = svaba_dir + '{sample}/{sample}.discordant.txt.gz',
        manta = manta_dir + '{sample}/results/variants/somaticSV.vcf.gz'
   output:
        'completed_files/{sample}.done'
   shell:
        'touch {output}'

rule rm_fq:
    input:
        gbam = bwa_dir + '{sample}-5.bam',
        gbai = bwa_dir + '{sample}-5.bam.bai',
        sbam = bwa_dir + '{sample}-1.bam',
        sbai = bwa_dir + '{sample}-1.bam.bai',
        donefile = 'completed_files/{sample}.done'
    output:
        'completed_files/{sample}_bams_removed'
    shell:
        'if [ -f {input.gbam} ]; then ' +
            'rm {input.gbam}; ' +
        'fi; ' + 
        'if [ -f {input.gbai} ]; then ' +
            'rm {input.gbai}; ' +
        'fi; ' + 
        'if [ -f {input.sbam} ]; then ' +
            'rm {input.sbam}; ' +
        'fi; ' + 
        'if [ -f {input.sbai} ]; then ' +
            'rm {input.sbai}; ' +
        'fi; ' + 
        'touch {output}'

