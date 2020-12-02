REF_PATH = '/Users/wckdouglas/Desktop/sars/ref'
PRIMER_SEQUENCE = REF_PATH + '/primers.fa'
PRIMER_COORDINATE = REF_PATH + '/primers.bed' 
PRIMER_FILE = REF_PATH + '/primers.txt' #see fgbio http://fulcrumgenomics.github.io/fgbio/tools/latest/TrimPrimers.html
BT_INDEX = PRIMER_SEQUENCE + '.1.bt2'

rule all:
    input:
        PRIMER_FILE, BT_INDEX

rule build_index:
    input:
        PRIMER_SEQUENCE

    output:
        BT_INDEX
    
    shell:
        'bowtie2-build {input} {input}'
        

rule make_primer:
    input:
        PRIMER = PRIMER_COORDINATE

    output:
        PRIMER = PRIMER_FILE

    run:
        import pandas as pd
        pd.read_csv(input.PRIMER,
                    usecols = [0,1,2,3], sep='\t',
                    names = ['chrom','start','end','name']) \
            .assign(coor = lambda d: d.start.astype(str) + ',' + d.end.astype(str))\
            .assign(strand = lambda d: d['name'].str.split('_',expand=True).iloc[:,-1]) \
            .assign(name = lambda d: d['name'].str.replace('_LEFT|_RIGHT',''))\
            .pipe(pd.pivot_table, 
                index=['chrom','name'], 
                columns = 'strand',
                values = 'coor', 
                aggfunc = lambda x: ','.join(x)) \
            .assign(left_start = lambda d: d.LEFT.str.split(',',expand=True).iloc[:,0].astype(int) + 1) \
            .assign(left_end = lambda d: d.LEFT.str.split(',',expand=True).iloc[:,1].astype(int) + 1) \
            .assign(right_start = lambda d: d.RIGHT.str.split(',',expand=True).iloc[:,0].astype(int) + 1) \
            .assign(right_end = lambda d: d.RIGHT.str.split(',',expand=True).iloc[:,1].astype(int) + 1) \
            .reset_index() \
            .filter(['chrom','left_start','left_end','right_start','right_end'])\
            .to_csv(output.PRIMER, sep='\t', index=False)

rule map:
    input:
        PRIMER_SEQUENCE

    params:
        REF = REF_PATH + '/Bavtpat1_complete.fa '

    output:
        PRIMER_COORDINATE

    shell:
        'bowtie2 -N1 -L 9 -x ${REF} -f {input} '\
        '| samtools view -b '\
        '| bedtools bamtobed -i - '\
        '> primers.bed'
