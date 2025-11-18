# MamPipe Docker Image for Google Cloud Batch
# Start with Miniconda3 base image
FROM continuumio/miniconda3:latest

# Avoid interactive prompts & speed up conda
ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_NO_PROGRESS_BARS=1
ENV MAMBA_ALLOW_CONDA_INSTALL=1
ENV PATH /opt/conda/envs/base_env/bin:$PATH

# Install Mamba and prepare environment
RUN conda install -y -n base -c conda-forge mamba && \
    mamba update -n base mamba && \
    mamba create -y -n base_env python=3.10 && \
    mamba install -y -n base_env \
        -c bioconda -c conda-forge \
        chopper \
        flye \
        medaka \
        busco \
        purge_dups \
        quast \
        gatk4 \
        ragtag \
        minimap2 \
        samtools \
        bwa \
        bcftools && \
    conda clean -afy

# Install NanoPlot via pip
RUN pip install NanoPlot

# Activate environment automatically
SHELL ["conda", "run", "-n", "base_env", "/bin/bash", "-c"]

# Set working directory
WORKDIR /workspace

# Set default command
CMD ["bash"]
