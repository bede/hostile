# Bioconda-based build against main branch
FROM condaforge/miniforge3:latest
RUN git clone https://github.com/bede/hostile.git
RUN sed -i 's/name: hostile/name: base/' hostile/environment.yml
RUN mamba env update -f hostile/environment.yml
RUN pip install ./hostile[dev]
RUN cd hostile && pytest

## Debian-based build against main branch
# FROM ubuntu:22.04
# RUN apt-get update && apt-get -y install bowtie2 minimap2 bedtools samtools python3-pip git
# RUN git clone https://github.com/bede/hostile.git
# RUN pip install ./hostile[dev]
# RUN hostile --version && bowtie2 --version && minimap2 --version && samtools --version && bedtools --version
# RUN cd hostile && pytest
