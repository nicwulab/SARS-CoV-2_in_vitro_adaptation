FROM ubuntu:18.04
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh;
RUN conda --version;

RUN conda config --set always_yes yes --set changeps1 no &&\
    conda info -a &&\
    conda list && \
    conda config --show-sources && \
    conda config --show && \
    conda config --add channels bioconda && \
    conda config --add channels default && \
    conda config --add channels anaconda && \
    conda config --add channels conda-forge;


RUN  conda install -c conda-forge mamba  
RUN mamba install -c bioconda -c anaconda python=3.6 \
        cutadapt bowtie2  \
        samtools snakemake pandas \
        varscan bamutil seqtk \
        mosdepth matplotlib \
        logomaker seaborn pytest && \
     pip install pysam; 


CMD [ "snakemake", "-s" ]
