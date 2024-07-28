bootstrap: docker
from: mambaorg/micromamba:1.5.8-bookworm-slim

%files 

    modules/local/R/environment.yml /opt

%post
    
    # install packages
    micromamba install -n base -f /opt/environment.yml
    micromamba clean -ay

    # need to get ps for Nextflow
    apt-get update
    apt-get install procps -y

%environment

    # use posix shell hooks
    eval "$(micromamba shell hook -s posix)"
    micromamba activate base

%labels
  Author yang.e@wehi.edu.au
  Version v0.1

%help
  This container is used in the SkewX pipeline developed by Quentin Gouil.
  It is built using the mambaor/micromamba container and installed packages
  can be viewed with micromamba list.
