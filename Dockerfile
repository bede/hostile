# Bioconda-based build against main branch
FROM condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive

COPY src/ hostile/
COPY environment.yml pyproject.toml README.md hostile/
WORKDIR /hostile

RUN sed -i 's/name: hostile/name: base/' environment.yml
RUN conda update conda
RUN conda env update -f environment.yml

RUN hostile --version && bowtie2 --version && minimap2 --version && samtools --version && bedtools --version
