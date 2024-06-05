bootstrap: docker
from: mambaorg/micromamba:1.5.8-bookworm-slim

%files 

    modules/local/R/environment.yml /opt

%post
    
    micromamba install -n base -f /opt/environment.yml -y
    micromamba install -n base -c conda-forge quarto
    micromamba clean -ay

%environment

    export PATH=/opt/conda/bin:$PATH
    export QUARTO_DENO=/opt/conda/bin/deno
