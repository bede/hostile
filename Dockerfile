FROM mambaorg/micromamba:2.0.2
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN cd /tmp && pip install '.[dev]'
RUN cd /tmp && pytest
