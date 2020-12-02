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

rule variant_calling:
    input:
        SORTED_BAM

    params:
        REF_FA = REF

    output:
        SNP_FILE
    
    shell:
        'samtools mpileup -f {params.REF_FA} {input} '\
        '| varscan pileup2snp '\
        ' --min-coverage 10 ' \
        ' --min-reads2 2 '\
        '--min-var-freq 0.01 '\
        '--min-freq-for-hom 0.75 '
        '--p-values 99e-02 ' \
        '> {output}'

rule sort:
    input:
        TRIMMED_BAM

    output:
        SORTED_BAM
    
    shell:
        'samtools sort {input} > {output};'\
        ' samtools index {output}'

rule trim_primer:
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


rule align:
    input:
        TRIMMED_FQ

    params:
        REF = REF

    output:
        BAM

    shell:
        'bowtie2 -x {params.REF} --interleaved {input} | samtools view -b  > {output}'

rule trim_adapt:
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
