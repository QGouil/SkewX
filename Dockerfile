FROM mambaorg/micromamba:1.5.8-bookworm-slim

COPY modules/local/R/environment.yml /opt

RUN micromamba install -n base -f /opt/environment.yml -y && micromamba clean -ay

# need to get ps for Nextflow, and install Quarto (conda's doesn't work for some reason)
RUN apt-get update && \
    apt-get install procps wget -y \
    wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.554/quarto-1.4.554-linux-amd64.deb -O /opt/quarto.deb \
    apt-get install /opt/quarto.deb -y \
    rm /opt/quarto.deb \
    apt-get clean

LABEL "Author"="yang.e@wehi.edu.au"
LABEL "Version"="v0.1"
