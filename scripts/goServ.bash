#!/bin/bash

#CODE1="cd pipelines/tbrc"
#CODE2="mkdir -p data/runs/"
#CODE3="rsync -a -P"

MUCID=" user@storage.example.org:/path/to/pipelines/tbrc/rawdata/"
VICTO=" user@storage.example.org:/path/to/pipelines/tbrc/rawdata/"

if [ $2 == storage2 ]; then CODE4=$MUCID; fi
if [ $2 == storage1 ]; then CODE4=$VICTO; fi

CONDA="TSV=data/runs/${1}/${4}/sample.tsv"

INPUTPATH="INTPT=data/runs/${1}/${4}/"

#CODE5="data/runs/"

#Anaconda path 
PATH="/path/to/scratch/anaconda3/condabin:/path/to/home/local/bin/:/usr/lib64/qt-3.3/bin:/path/to/home/perl5/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/path/to/home/bin/"

ssh dtn.hpc.example.org "cd pipelines/tbrc; mkdir -p data/runs/${1}; rsync -a -P${CODE4}${3} data/runs/${1}"
rsync -a -P runs/${1}/${4}/ dtn.hpc.example.org:pipelines/tbrc/data/runs/${1}/${4}/

ssh user@login.hpc.example.org "cd pipelines/tbrc;  export PATH=${PATH}; echo $PATH; ${CONDA}; ${INPUTPATH}; export TSV; export INTPT; conda activate snakemake; ./snakemake.sh"
