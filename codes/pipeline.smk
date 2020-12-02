# only variable needed to change
PROJECT_PATH='/Users/wckdouglas/Desktop/sars'


REF = PROJECT_PATH + '/ref/Bavtpat1_complete.fa'
PRIMERS = PROJECT_PATH + '/ref/primers.txt'
R1 = PROJECT_PATH + '/data/{SAMPLENAME}_L001_R1_001.fastq.gz'
R2 = PROJECT_PATH + '/data/{SAMPLENAME}_L001_R2_001.fastq.gz'
SAMPLENAMES, = glob_wildcards(R1)
RESULT_PATH = PROJECT_PATH + '/results/{SAMPLENAME}'
TRIMMED_FQ = RESULT_PATH + '/trimmed.fq.gz'
MERGED_FQ = RESULT_PATH + '/merged.fq.gz'
BAM = RESULT_PATH + '/aligned.bam'
TRIMMED_BAM = RESULT_PATH + '/trimmed.bam'
SORTED_BAM = RESULT_PATH + '/sorted.bam'
SNP_FILE = RESULT_PATH + '/variants.snp'

rule all:
    input:
        expand(SNP_FILE, SAMPLENAME = SAMPLENAMES)

rule variant_calling_with_varscan:
    input:
        SORTED_BAM

    params:
        REF_FA = REF

    output:
        SNP_FILE
    
    shell:
        'samtools mpileup --fasta-ref {params.REF_FA} '\
        '--max-depth 50000 --min-MQ 30 --min-BQ 30  {input} '\
        '| varscan pileup2cns '\
        ' --min-coverage 10 ' \
        ' --min-reads2 2 '\
        '--min-var-freq 0.01 '\
        '--min-freq-for-hom 0.75 '
        '--p-value 0.05 ' \
        '> {output}'

rule sort_bam_with_samtools:
    input:
        TRIMMED_BAM

    output:
        SORTED_BAM
    
    shell:
        'samtools sort {input} > {output};'\
        ' samtools index {output}'

rule trim_primers_from_alignment_with_fgbio:
    input:
        BAM

    params:
        REF_FA = REF,
        PRIMERS = PRIMERS

    output:
        TRIMMED_BAM

    shell:
        'fgbio TrimPrimers '\
        '-i {input} ' \
        '-o {output} '\
        '-p {params.PRIMERS} '\
        '-r {params.REF_FA} '


rule align_with_bowtie:
    input:
        TRIMMED_FQ

    params:
        REF = REF

    output:
        BAM

    shell:
        'bowtie2 -x {params.REF} --interleaved {input} | samtools view -b  > {output}'

rule trim_adapter_with_cutadapt:
    input:
        FQ1 = R1,
        FQ2 = R2
        
    output:
        TRIMMED_FQ

    shell:
        'seqtk mergepe {input.FQ1} {input.FQ2} '\
        '| cutadapt -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG '\
        '-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC '\
        '--interleaved --minimum-length 50 '\
        '-o {output} -'
