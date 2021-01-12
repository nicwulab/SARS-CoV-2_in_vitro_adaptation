# only variable needed to change
PROJECT_PATH='/Users/wckdouglas/Desktop/sars'

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)

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
DEPTH_FILE = RESULT_PATH + '/coverage.per-base.bed.gz'
SNP_FILE = RESULT_PATH + '/variants.snp'
SEQ_LOGO = RESULT_PATH + '/MBCS_seqlogo.png'
FREQ_FILE = RESULT_PATH + '/MBCS_freq.tsv'
COVERAGE_PNG = PROJECT_PATH + '/results/coverage.png'
COVERAGE_STAT = PROJECT_PATH + '/results/coverage_stat.csv'

rule all:
    input:
        expand(SNP_FILE, SAMPLENAME = SAMPLENAMES),
        expand(SEQ_LOGO, SAMPLENAME = SAMPLENAMES),
        expand(FREQ_FILE, SAMPLENAME = SAMPLENAMES),
        COVERAGE_PNG, COVERAGE_STAT

rule plotDepth:
    input:
        expand(DEPTH_FILE, SAMPLENAME = SAMPLENAMES)

    output:
       FIG =  COVERAGE_PNG,
       TAB = COVERAGE_STAT

    run:
        def read_bed(bed):
            return pd.read_csv(bed, names = ['chrom','start','end','read_coverage'], sep='\t') \
                .assign(samplename = os.path.basename(bed).split('.')[0]) \
                .assign(log_cov = lambda d: d.read_coverage.transform(np.log))
        
        
        dfs = map(read_bed, input)
        df = pd.concat(dfs)
        
        #make table
        df\
            .groupby('samplename', as_index=False)\
            .agg({'read_coverage':'median'})\
            .rename(columns = {'read_coverage':'median read coverage'})\
            .merge(df\
                    .groupby('samplename', as_index=False)\
                    .agg({'read_coverage':'mean'})\
                    .rename(columns = {'read_coverage':'mean read coverage'}))\
            .to_csv(output.TAB, index=False)
    
        #make figure
        p = sns.FacetGrid(data = df, col_wrap = 5, col = 'samplename')
        p.map(plt.plot, 'end', 'log_cov')
        p.set_titles(col_template = '{col_name}')
        sns.despine()
        p.set(xlabel = 'Read coverage (log)', ylabel = 'Genomic position')
        fig.savefig(output.FIG, bbox_inches='tight')




rule cal_depth:
    input:
        SORTED_BAM

    params:
        PREFIX = RESULT_PATH + '/coverage'

    output:
        DEPTH_FILE 
    
    shell:
        'mosdepth {params.PREFIX} {input}'
        

rule cal_frequency:
    input:
        SORTED_BAM

    output:
        SEQ_LOGO = SEQ_LOGO, 
        FREQ_FILE = FREQ_FILE

    shell:
        'python extract_MBCS.py {input} {output.FREQ_FILE} {output.SEQ_LOGO}'


rule variant_calling_with_varscan:
    input:
        SORTED_BAM

    params:
        REF_FA = REF

    output:
        SNP_FILE
    
    shell:
        'samtools mpileup --excl-flags 2048 --excl-flags 256  '\
        '--fasta-ref {params.REF_FA} '\
        '--max-depth 50000 --min-MQ 30 --min-BQ 30  {input} '\
        '| varscan pileup2cns '\
        ' --min-coverage 10 ' \
        ' --min-reads2 2 '\
        '--min-var-freq 0.01 '\
        '--min-freq-for-hom 0.75 '
        '--p-value 0.05 --variants 1 ' \
        '> {output}'

rule sort_bam_with_samtools:
    input:
        TRIMMED_BAM

    output:
        SORTED_BAM
    
    shell:
        'samtools sort {input} > {output};'\
        ' samtools index {output}'

rule trim_primers_from_alignment_with_bamutils:
    input:
        BAM

    params:
        REF_FA = REF,
        PRIMERS = PRIMERS

    output:
        TRIMMED_BAM

    shell:
        'bam trimbam {input} - -L 30 -R 0 --clip '\
        '| samtools fixmate - - '\
        '| samtools calmd -Q - {params.REF_FA} '\
        '> {output} '


rule align_with_bowtie:
    input:
        TRIMMED_FQ

    params:
        REF = REF

    output:
        BAM

    shell:
        'bowtie2 -x {params.REF} '\
        '--no-discordant --dovetail --no-mixed --maxins 2000 ' \
        '--interleaved {input} --mm '\
        '| samtools view -bF 4  > {output}'

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
