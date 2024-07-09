# Bioconda-based build against main branch
FROM condaforge/miniforge3:latest
RUN git clone https://github.com/EIT-Pathogena/hostile-eit.git
RUN mamba env update -n base -f hostile/environment.yml
RUN pip install ./hostile[dev]
RUN cd hostile && pytest
