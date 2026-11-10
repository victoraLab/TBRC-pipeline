# TBRC images. Two targets:
#
#   --target backend   pipeline runtime only: snakemake + every tool + R + clonality
#   --target unified   backend plus the Shiny app
#
# `./tbrc setup` builds what is missing. By hand, from the repository root:
#
#   docker build --target backend -t tbrc-backend:arm64 deploy/docker   # `./tbrc run --local` looks for this name
#   docker build --target unified -t tbrc-unified:arm64 deploy/docker   # the web app (docker-compose.yml) runs this
#
# The first build takes about 15 minutes and needs internet (conda channels and GitHub); the
# second reuses the layers of the first.
#
# Tested on one Apple Silicon Mac; not tested on Linux or Intel. The tags say arm64 because that
# is what was built and run (Docker Desktop, linux/arm64). Nothing here selects a platform, so
# Docker builds for the machine it runs on, which may not work for another architecture.
#
# Neither the workflow nor the app is baked in. Both are bind-mounted from the clone, so code
# updates arrive with `git pull` and never need an image rebuild. Rebuild only to change a
# tool version.
#
# The nanopore tools (abpoa, racon, minimap2, mafft, spoa, seqkit) are installed here and sit
# on PATH. backend/config/pipelines.yaml's nanopore tools.env_bin names a conda environment on
# a cluster; that path does not exist here and is prepended to PATH harmlessly.
#
# Versions are pinned to what the lab's cluster runs (checked 2026-11-06), so a run here and a
# run there execute the same tools. Unpinned packages float, so a rebuild months from now may
# resolve other versions than the ones the pins were checked against.

# ---------------------------------------------------------------- base
FROM condaforge/miniforge3:25.3.0-3 AS base
ENV DEBIAN_FRONTEND=noninteractive

# gawk and perl: fastx_barcode_splitter.pl and split_fasta.sh; rsync/ssh: only for a
# cluster or storage server; jq parses server_config.json; libltdl7 is a pandaseq runtime
# dependency its conda package does not pull in.
RUN apt-get update && apt-get install -y --no-install-recommends \
        zip unzip rsync openssh-client jq procps ca-certificates gawk perl libltdl7 \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------- backend
FROM base AS backend

# blast is pinned because 2.17.0 has no linux-aarch64 build; unpinned, the solver holds
# at 2.16.0 here and a rebuild elsewhere would diverge. cutadapt and pandas are pinned to the
# cluster's versions: unpinned they resolved to cutadapt 5.2 and pandas 3.1 (a major version
# with behaviour changes) on 2026-11-06. igblast resolves to 1.22.0, the cluster's version.
RUN mamba install -y -n base -c conda-forge -c bioconda \
        "snakemake-minimal=7.32.4" \
        "pandaseq=2.11" \
        "fastx_toolkit=0.0.14" \
        "cutadapt=4.8" \
        "igblast=1.22.0" \
        "vsearch=2.31.0" \
        "blast=2.16.0" \
        "pandas=2.1.2" \
    && mamba clean -afy

# Nanopore pipeline (tiago_nanopore): same versions as the lab's tbrc_nanopore env.
RUN mamba install -y -n base -c conda-forge -c bioconda \
        "abpoa=1.5.6" "racon=1.5.0" "minimap2=2.31" "mafft=7.526" "spoa=4.1.5" "seqkit=2.13.0" \
    && mamba clean -afy

# QC, plate figures and clonality (clonality is the lab's own R package, public on GitHub).
RUN mamba install -y -n base -c conda-forge -c bioconda \
        r-base r-ggplot2 r-dplyr r-readr r-remotes \
    && Rscript -e "options(Ncpus=4); \
        install.packages('ggplate', repos='https://cloud.r-project.org'); \
        remotes::install_github('victoraLab/clonality', upgrade='never')" \
    && Rscript -e "stopifnot(requireNamespace('clonality'), requireNamespace('ggplate'))" \
    && mamba clean -afy

WORKDIR /work
CMD ["snakemake", "--cores", "4", "--printshellcmds"]

# ---------------------------------------------------------------- unified
# The app plus the pipeline. The app's capability check finds the tools on PATH.
FROM backend AS unified
RUN mamba install -y -n base -c conda-forge \
        r-shiny r-shinythemes r-shinyalert r-spscomps \
        r-dt r-jsonlite r-shinywidgets \
    && mamba clean -afy
RUN mkdir -p /root/.ssh && chmod 700 /root/.ssh
WORKDIR /srv/tbrc
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/srv/tbrc', host='0.0.0.0', port=3838, launch.browser=FALSE)"]
