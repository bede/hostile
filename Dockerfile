FROM ubuntu:22.04
RUN apt-get update && apt-get -y install bowtie2 minimap2 bedtools samtools python3-pip git
RUN git clone https://github.com/bede/hostile.git
RUN pip install ./hostile[dev]
RUN hostile --version && bowtie2 --version && minimap2 --version && samtools --version && bedtools --version
RUN cd hostile && pytest


## Hermetic bioconda build
# FROM continuumio/miniconda3:latest
# WORKDIR /hostile
# COPY . /hostile
# RUN conda env update -f environment.yml
# RUN pip install ./
# RUN hostile --version && bowtie2 --version && minimap2 --version && samtools --version && bedtools --version


# FROM condaforge/miniforge3:latest
# RUN mamba install -y -c conda-forge -c bioconda minimap2


# FROM continuumio/miniconda3:latest
# RUN conda install -y -c conda-forge -c bioconda minimap2
