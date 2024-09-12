# Bioconda-based build against main branch
FROM condaforge/miniforge3:24.1.2-0
RUN git clone https://github.com/bede/hostile.git
RUN mamba env update -n base -f hostile/environment.yml
RUN pip install ./hostile[dev]
RUN cd hostile && pytest
