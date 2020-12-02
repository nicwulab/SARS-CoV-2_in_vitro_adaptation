# Dependencies #

* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org/)
* [fgbio](https://github.com/fulcrumgenomics/fgbio)
* [varscan](http://dkoboldt.github.io/varscan/)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [pandas](https://pandas.pydata.org/)

## Installation ##

Easiest way to install everything is through [miniconda](https://docs.conda.io/en/latest/miniconda.html)

And the following to install the needed packages:
```
conda install -c bioconda -c anaconda python=3.6 cutadapt bowtie2 samtools fgbio snakemake pandas
```

## Steps ##

1. build reference:
    - go into  ```/code/``` folder
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
| BavtPat1.fa |    2,619 | T   | Y    |     98 |      2 |    2.00 |        1 |        1 |    40 |    34 |   0.98 |     True |     True |         98 |           0 |          2 |           0 | C         |
| BavtPat1.fa |    3,688 | C   | Y    |    795 |     79 |    9.04 |        2 |        1 |    70 |    72 |   0.98 |     True |     True |        794 |           1 |         79 |           0 | T         |
| BavtPat1.fa |    3,813 | A   | M    |  1,003 |     22 |    2.15 |        2 |        1 |    38 |    37 |   0.98 |     True |     True |        325 |         678 |          0 |          22 | C         |
| BavtPat1.fa |    4,603 | T   | Y    |     88 |      2 |    2.22 |        2 |        2 |    55 |    24 |   0.98 |     True |     True |         51 |          37 |          1 |           1 | C         |
| BavtPat1.fa |    4,719 | C   | Y    |     41 |      2 |    4.55 |        2 |        1 |    70 |    51 |   0.98 |     True |     True |         40 |           1 |          2 |           0 | T         |
| BavtPat1.fa |    4,720 | G   | K    |     41 |      2 |    4.55 |        1 |        1 |    72 |    60 |   0.98 |     True |     True |         41 |           0 |          2 |           0 | T         |
| BavtPat1.fa |    4,728 | G   | S    |     45 |      2 |    4.17 |        2 |        1 |    72 |    40 |   0.98 |     True |     True |         44 |           1 |          2 |           0 | C         |
| BavtPat1.fa |    4,730 | T   | W    |     45 |      2 |    3.92 |        2 |        1 |    72 |    62 |   0.98 |     True |     True |         43 |           2 |          2 |           0 | +C        |
```






