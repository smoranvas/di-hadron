#!/bin/bash

A=$1
i=$2
export RIVET_CONTAINER_NAME=${PWD}/pythia-eic-tutorial_latest.sif
alias rivet='apptainer run $RIVET_CONTAINER_NAME rivet'
cp *.so /volatile/clas12/spaul/beagle_output/;
cd /volatile/clas12/spaul/beagle_output/; apptainer run $RIVET_CONTAINER_NAME rivet -a pion_proton_correlations --pwd -o e${A}_${\
i}.yoda e${A}_${i}.hepmc;
