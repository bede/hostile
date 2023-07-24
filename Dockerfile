FROM ubuntu:22.04
RUN apt-get update && apt-get -y install bowtie2 minimap2 bedtools samtools python3-pip git
RUN pip install hostile[dev]
RUN hostile --version && bowtie2 --version && minimap2 --version && samtools --version && bedtools --version
RUN git clone https://github.com/bede/hostile.git && cd hostile && pytest
