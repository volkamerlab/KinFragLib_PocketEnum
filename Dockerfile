FROM mambaorg/micromamba:1.4.3
USER root
COPY --chown=kabu00002:kabu00002 environment.yml /tmp/env.yaml
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN apt-get update && \
        apt-get -y install sudo \
        pip 
RUN apt-get install -y git
RUN sudo apt-get install -y libglib2.0-0

# Create a working directory
RUN mkdir /main
RUN mkdir /main/home
WORKDIR /main/home

RUN mkdir /main/home/.cache
RUN mkdir /main/home/.local
RUN mkdir /main/home/.local/share
RUN mkdir /main/home/.cache/wandb
RUN mkdir /main/home/.local/share/wandb

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

RUN git clone https://github.com/volkamerlab/KinFragLib.git && \
        pip install -e KinFragLib
