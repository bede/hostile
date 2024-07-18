# Bioconda-based build against main branch
FROM condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive

COPY src/ hostile/
COPY environment.yml pyproject.toml README.md hostile/
WORKDIR /hostile

RUN conda update conda
RUN conda env update --name base -f environment.yml

RUN hostile --version && bowtie2 --version && minimap2 --version && samtools --version && bedtools --version
