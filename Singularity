bootstrap: docker
from: mambaorg/micromamba:1.5.8-bookworm-slim

%files 

    modules/local/R/environment.yml /opt

%post
    
    micromamba install -n base -f /opt/environment.yml -y
    micromamba clean -ay

    # need to get ps for Nextflow, and install Quarto (conda's doesn't work for some reason)
    apt-get update
    apt-get install procps wget -y
    wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.554/quarto-1.4.554-linux-amd64.deb -O /opt/quarto.deb
    apt-get install /opt/quarto.deb -y
    rm /opt/quarto.deb

%environment

    export PATH=/opt/conda/bin:$PATH

%labels
  Author yang.e@wehi.edu.au
  Version v0.1

%help
  This container is used in the SkewX pipeline developed by Quentin Gouil.
  It is built using the mambaor/micromamba container and installed packages
  can be viewed with micromamba list.
