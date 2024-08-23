FROM mambaorg/micromamba:1.5.8-noble

LABEL author="Jinlong Ru"

USER root

# Activate conda env during docker build
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Set PATH manually, because nextflow doesn't activate base env by default.
ENV PATH=/opt/conda/bin:$PATH

# Copy conda env files
COPY ./docker/tpd-base/envs*.yml /tmp

# Install packages
RUN --mount=type=cache,target=/opt/conda/pkgs micromamba install -n base -f /tmp/envs-base.yml -y

# clean
RUN micromamba clean --all --yes