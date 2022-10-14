#!/bin/bash

#CODE1="rsync -a -P"

FOLDER=runs/${1}/${4%.*}

mkdir -p ${FOLDER}
cp ${2}  ${FOLDER}/${4}
cp ${3}  ${FOLDER}/${5}

#${CODE1} runs/${1} dtn.hpc.example.org:pipelines/tbrc/data/runs/
rsync -a -P runs/${1} dtn.hpc.example.org:pipelines/tbrc/data/runs/


