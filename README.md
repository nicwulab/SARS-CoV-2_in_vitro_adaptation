# SARS-CoV-2 cell culture-adaptive mutations #

This is the code repository for the manuscript [Human airway cells prevent SARS-CoV-2 multibasic cleavage site cell culture adaptation
](https://www.biorxiv.org/content/10.1101/2021.01.22.427802v1).

# Dependencies #

* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org/)
* [varscan](http://dkoboldt.github.io/varscan/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [bamutil](https://genome.sph.umich.edu/wiki/BamUtil:_trimBam)
* [mosdepth](https://github.com/brentp/mosdepth)
* [pandas](https://pandas.pydata.org/)
* [matplotlib](https://matplotlib.org/)
* [seaborn](https://seaborn.pydata.org/)

## Installation ##

Easiest way to install everything is through [miniconda](https://docs.conda.io/en/latest/miniconda.html)

And the following to install the needed packages:
```
conda create -n SARS -c bioconda -c anaconda python=3.6 cutadapt  \
            bowtie2 samtools snakemake pandas varscan bamutil seqtk   \
            mosdepth matplotlib pandas pysam logomaker seaborn
```

Before running the analysis, do:
```
conda activate SARS
```

This requires doing ```conda init``` first.

## Steps ##

1. build reference:
    - go into  ```/code/``` folder
    - set the ```REF_PATH``` variable in ```make_ref.smk``` file accordingly
    - do ```snakemake -s make_ref.smk```
2. Download all fastq files into ```/data/``` folder
3. Run analysis:
    - go into ```/code/``` folder
    - set the ```PROJECT_PATH``` variable in ```pipeline.smk``` file accordingly
    - do: ```snakemake -s pipeline.smk``` to excute the analysis
4. For each sample, you will get the following files:
    ```
    results
    ├── P1_S1
    │   ├── aligned.bam
    │   ├── sorted.bam
    │   ├── sorted.bam.bai
    │   ├── trimmed.bam
    │   ├── trimmed.fq.gz
    │   └── variants.snp
    ```

The excuted workflow is as following: 

![workflow](https://github.com/wckdouglas/SARS_CoV2_mutation/blob/main/codes/pipeline.png?raw=true)

## Result SNP file (variants.snp) ##

```
| Chrom       | Position | Ref | Cons | Reads1 | Reads2 | VarFreq | Strands1 | Strands2 | Qual1 | Qual2 | Pvalue | MapQual1 | MapQual2 | Reads1Plus | Reads1Minus | Reads2Plus | Reads2Minus | VarAllele |
| ----------- | -------- | --- | ---- | ------ | ------ | ------- | -------- | -------- | ----- | ----- | ------ | -------- | -------- | ---------- | ----------- | ---------- | ----------- | --------- |
| BavtPat1.fa |    3,688 | C   | Y    |    795 |     79 |    9.04 |        2 |        1 |    70 |    72 | 0.000… |     True |     True |        794 |           1 |         79 |           0 | T         |
| BavtPat1.fa |    3,813 | A   | M    |  1,003 |     22 |    2.15 |        2 |        1 |    38 |    37 | 0.000… |     True |     True |        325 |         678 |          0 |          22 | C         |
| BavtPat1.fa |    5,854 | T   | K    |  7,363 |    122 |    1.62 |        2 |        1 |    35 |    21 | 0.000… |     True |     True |      7,334 |          29 |        122 |           0 | G         |
| BavtPat1.fa |   10,086 | G   | R    |  7,289 |     87 |    1.17 |        2 |        1 |    35 |    21 | 0.000… |     True |     True |      7,199 |          90 |         87 |           0 | A         |
| BavtPat1.fa |   14,120 | C   | Y    |    101 |      8 |    7.34 |        2 |        2 |    63 |    64 | 0.003… |     True |     True |         82 |          19 |          6 |           2 | T         |
| BavtPat1.fa |   18,228 | T   | Y    |    238 |     58 |   19.59 |        2 |        1 |    41 |    19 | 0.000… |     True |     True |        214 |          24 |         58 |           0 | C         |
```


The columns represents:
```
OUTPUT
	Tab-delimited SNP calls with the following columns:
	Chrom		chromosome name
	Position	position (1-based)
	Ref		reference allele at this position
	Cons		Consensus genotype of sample in IUPAC format.
	Reads1		reads supporting reference allele
	Reads2		reads supporting variant allele
	VarFreq		frequency of variant allele by read count
	Strands1	strands on which reference allele was observed
	Strands2	strands on which variant allele was observed
	Qual1		average base quality of reference-supporting read bases
	Qual2		average base quality of variant-supporting read bases
	Pvalue		Significance of variant read count vs. expected baseline error
	MapQual1	Average map quality of ref reads (only useful if in pileup)
	MapQual2	Average map quality of var reads (only useful if in pileup)
	Reads1Plus	Number of reference-supporting reads on + strand
	Reads1Minus	Number of reference-supporting reads on - strand
	Reads2Plus	Number of variant-supporting reads on + strand
	Reads2Minus	Number of variant-supporting reads on - strand
	VarAllele	Most frequent non-reference allele observed 
```




