FROM continuumio/miniconda3

RUN conda config --set always_yes yes --set changeps1 no;\
    conda info -a;\
    conda list;\
    conda config --show-sources;\
    conda config --show;

RUN  conda install -c conda-forge mamba; \
     mamba install -c bioconda -c anaconda python=3.6 \
        cutadapt bowtie2  \
        samtools snakemake pandas \
        varscan bamutil seqtk \
        mosdepth matplotlib pandas \ 
        logomaker seaborn pytest; \
     pip install pysam; 

COPY . /opt/pipeline

CMD [ "snakemake", "-s" ,"/opt/pipeline/codes/pipeline.smk"]


