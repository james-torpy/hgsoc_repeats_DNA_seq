# Run command:
# snakemake --reason --cores 150 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N cnv_ident.smk -wd '/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/DNA-seq/logs' -b y -j y -V' -j 10

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script remaps bams to hg38 and uses SvABA and Manta to identify breakpoints in HGSOC genomes ###

# define directories:
project_name = 'hgsoc_repeats'
exp_name = 'DNA-seq'

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/' + exp_name + '/'
results_dir = project_dir + 'results/'
genome_dir = project_dir + 'genome/'
script_dir = project_dir + 'scripts/'

bam_dir = 'raw_files/bams/'
temp_dir = project_dir + bam_dir + '/temp/'
fq_dir = 'raw_files/fastq/'
bwa_dir = 'results/bwa/'
svaba_dir = 'results/svaba/'
manta_dir = 'results/manta/'
manta_bin = '/g/data1a/ku3/jt3341/local/lib/manta-1.5.0/bin/'

SAMPLES = list([
#    "AOCS-063-sub"
    "AOCS-063", "AOCS-064", "AOCS-065", "AOCS-075", "AOCS-076" 
#    , "AOCS-077", "AOCS-078", "AOCS-080", 
#    "AOCS-083", "AOCS-084", "AOCS-085", "AOCS-086", 
#    "AOCS-090", "AOCS-091", "AOCS-092", "AOCS-093", 
#    "AOCS-094", "AOCS-095", "AOCS-107", "AOCS-108", 
#    "AOCS-109", "AOCS-111", "AOCS-112", "AOCS-113", 
#    "AOCS-114", "AOCS-115", "AOCS-116", "AOCS-122", 
#    "AOCS-123", "AOCS-124", "AOCS-125", "AOCS-126", 
#    "AOCS-128", "AOCS-130", "AOCS-131", "AOCS-133", 
#    "AOCS-137", "AOCS-143", "AOCS-144", "AOCS-145", 
#    "AOCS-146", "AOCS-147", "AOCS-148", "AOCS-149", 
#    "AOCS-152", "AOCS-153"
])

rule all:
    input:
        expand(
            'results/star/GC/{sample}/bams_removed',
            sample=SAMPLES
        )

rule transfer:
    output:
        fq1='raw_files/bams/{sample}-5.bam.gz',
        fq2='raw_files/bams/{sample}-1.bam.gz'
    threads: 2
    shell:
        'mkdir -p logs/transfer; ' + 
        'cd logs/transfer; ' + 
        ' ../../scripts/1.transfer_bam.sh' + 
        ' {wildcards.sample}' +
            ' 2> {wildcards.sample}.transfer.errors'

rule gunzip1:
    input:
        bam_dir + '{sample}-5.bam.gz'
    output:
        bam_dir + '{sample}-5.bam'
    threads: 8
    shell:
        'module load phuluu/pigz/2.3.4; ' +
        'mkdir -p logs/gunzip1; ' + 
        'cd logs/gunzip1; ' + 
        'pigz -d ../../{input}' +
            ' 2> {wildcards.sample}.gunzip1.errors'

rule gunzip2:
    input:
        bam_dir + '{sample}-1.bam.gz'
    output:
        bam_dir + '{sample}-1.bam'
    threads: 8
    shell:
        'module load phuluu/pigz/2.3.4; ' +
        'mkdir -p logs/gunzip2; ' + 
        'cd logs/gunzip2; ' + 
        'pigz -d ../../{input}' +
            ' 2> {wildcards.sample}.gunzip2.errors'

rule sort1a:
    input:
        bam_dir + '{sample}-5.bam'
    output:
        bam_dir + '{sample}-5.sorted.by.name.bam'
    threads: 10
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/sort1a; ' + 
        'cd logs/sort1a; ' + 
        'samtools sort -n -@ 6 ../../{input} -o ../../{output}' +
            ' 2> {wildcards.sample}.sort1a.errors; ' +
        'if [ -f {output} ]; then ' +
            'rm ../../{input}; ' +
        'fi'
        

rule sort1b:
    input:
        bam_dir + '{sample}-1.bam'
    output:
        bam_dir + '{sample}-1.sorted.by.name.bam'
    threads: 10
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/sort1b; ' + 
        'cd logs/sort1b; ' + 
        'samtools sort -n -@ 6 ../../{input} -o ../../{output}' +
            ' 2> {wildcards.sample}.sort1b.errors; ' +
        'if [ -f {output} ]; then ' +
            'rm ../../{input}; ' +
        'fi'

rule bam2fastq1:
    input:
        bam_dir + '{sample}-5.sorted.by.name.bam'
    output:
        fq1 = fq_dir + '{sample}-5-R1.fq',
        fq2 = fq_dir + '{sample}-5-R2.fq',
        fq_single = fq_dir + '{sample}-5-singles.fq'
    threads: 8
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/bam2fastq1; ' + 
        'cd logs/bam2fastq1; ' + 
        'samtools bam2fq -1 ../../{output.fq1} ' +
        ' -2 ../../{output.fq2} -@ 15 -s ../../{output.fq_single} ../../{input}' + 
                ' 2> {wildcards.sample}.bam2fastq1.errors; ' + 
        'if [ -f {output.fq1} ]; then ' +
            'if [ -f {output.fq2} ]; then ' +
                'rm ../../{input}; ' +
            'fi; ' +
        'fi'

rule bam2fastq2:
    input:
        bam_dir + '{sample}-1.sorted.by.name.bam'
    output:
        fq1 = fq_dir + '{sample}-1-R1.fq',
        fq2 = fq_dir + '{sample}-1-R2.fq',
        fq_single = fq_dir + '{sample}-1-singles.fq'
    threads: 8
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/bam2fastq2; ' + 
        'cd logs/bam2fastq2; ' + 
        'samtools bam2fq -1 ../../{output.fq1} ' +
        ' -2 ../../{output.fq2} -@ 15 -s ../../{output.fq_single} ../../{input}' + 
            ' 2> {wildcards.sample}.bam2fastq2.errors; ' +
        'if [ -f {output.fq1} ]; then ' +
            'if [ -f {output.fq2} ]; then ' +
                'rm ../../{input}; ' +
            'fi; ' +
        'fi'

rule bwa_align1:
    input:
        fq1 = fq_dir + '{sample}-5-R1.fq',
        fq2 = fq_dir + '{sample}-5-R2.fq',
        fq_single = fq_dir + '{sample}-5-singles.fq',
        ind = genome_dir + 'GRCh38.primary_assembly.genome.fa.amb'
    output:
        bwa_dir + '{sample}-5.sam'
    threads: 15
    shell:
        'module load marsmi/bwa/0.7.17; ' + 
        'mkdir -p logs/bwa_align1; ' + 
        'cd logs/bwa_align1; ' + 
        'bwa mem -t 10 ' + genome_dir + 
            'GRCh38.primary_assembly.genome.fa ../../{input.fq1} ../../{input.fq2}' + 
            '> ../../{output}; ' +
            'if [ -f {output} ]; then ' + 
               	'rm ../../{input.fq1}; ' + 
               	'rm ../../{input.fq2}; ' + 
               	'rm ../../{input.fq_single}; ' + 
            'fi'

rule bwa_align2:
    input:
        fq1 = fq_dir + '{sample}-1-R1.fq',
        fq2 = fq_dir + '{sample}-1-R2.fq',
        fq_single = fq_dir + '{sample}-1-singles.fq',
        ind = genome_dir + 'GRCh38.primary_assembly.genome.fa'
    output:
        bwa_dir + '{sample}-1.sam'
    threads: 15
    shell:
        'module load marsmi/bwa/0.7.17; ' + 
        'mkdir -p logs/bwa_align2; ' + 
        'cd logs/bwa_align2; ' + 
        'bwa mem -t 10 ' + genome_dir + 
            'GRCh38.primary_assembly.genome.fa ../../{input.fq1} ../../{input.fq2}' + 
            '> ../../{output}; ' +
            'if [ -f {output} ]; then ' + 
            	'rm ../../{input.fq1}; ' + 
            	'rm ../../{input.fq2}; ' + 
            	'rm ../../{input.fq_single}; ' + 
            'fi'

rule bam1:
    input:
        bwa_dir + '{sample}-5.sam'
    output:
        bwa_dir + '{sample}-5.unsorted.bam'
    threads: 15
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/bam1; ' + 
        'cd logs/bam1; ' + 
        'samtools view -bh  -@ 15 ../../{input} > ../../{output}; ' + 
        'if [ -f {output} ]; then ' + 
            'rm ../../{input}; ' + 
        'fi'

rule bam2:
    input:
        bwa_dir + '{sample}-1.sam'
    output:
        bwa_dir + '{sample}-1.unsorted.bam'
    threads: 15
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/bam2; ' + 
        'cd logs/bam2; ' + 
        'samtools view -bh -@ 15 ../../{input} > ../../{output}; ' + 
        'if [ -f {output} ]; then ' + 
            'rm ../../{input}; ' +
        'fi'

rule sort2a:
    input:
        bwa_dir + '{sample}-5.unsorted.bam'
    output:
        bwa_dir + '{sample}-5.bam'
    threads: 10
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/sort2a; ' + 
        'cd logs/sort2a; ' + 
        'samtools sort -@ 10 ../../{input} -o ../../{output}' +
            ' 2> {wildcards.sample}.sort2a.errors; ' +
        'if [ -f {output} ]; then ' + 
            'rm ../../{input}; ' +
        'fi'

rule sort2b:
    input:
        bwa_dir + '{sample}-1.unsorted.bam'
    output:
        bwa_dir + '{sample}-1.bam'
    threads: 10
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/sort2b; ' + 
        'cd logs/sort2b; ' + 
        'samtools sort -@ 10 ../../{input} -o ../../{output}' +
            ' 2> {wildcards.sample}.sort2b.errors; '
        'if [ -f {output} ]; then ' + 
            'rm ../../{input}; ' +
        'fi'

rule index1:
    input:
        bwa_dir + '{sample}-5.bam'
    output:
        bwa_dir + '{sample}-5.bam.bai'
    threads: 15
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/index1; ' + 
        'cd logs/index1; ' + 
        'samtools index -@ 15 ../../{input}' +
            ' 2> {wildcards.sample}.index1.errors'

rule index2:
    input:
        bwa_dir + '{sample}-1.bam'
    output:
        bwa_dir + '{sample}-1.bam.bai'
    threads: 15
    shell:
        'module load briglo/samtools/1.9; ' +
        'mkdir -p logs/index2; ' + 
        'cd logs/index2; ' + 
        'samtools index -@ 15 ../../{input}' +
            ' 2> {wildcards.sample}.index2.errors'

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
        'svaba run -t {input.sbam} -n {input.gbam} -G ' + 
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
        bam_dir + '{sample}.done'
   shell:
        'touch {output}'

rule rm_fq:
    input:
        gbam = bwa_dir + '{sample}-5.bam',
        gbai = bwa_dir + '{sample}-5.bam.bai',
        sbam = bwa_dir + '{sample}-1.bam',
        sbai = bwa_dir + '{sample}-1.bam.bai'
    output:
        'raw_files/bams/{sample}_bams_removed'
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

