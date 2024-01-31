#!/bin/bash

export RIVET_CONTAINER_NAME=${PWD}/pythia-eic-tutorial_latest.sif
alias rivet='apptainer $RIVET_CONTAINER_NAME rivet'
alias rivet-mkanalysis='apptainer run $RIVET_CONTAINER_NAME rivet-mkanalysis'
alias rivet-build='apptainer run $RIVET_CONTAINER_NAME rivet-build'
alias rivet-mkhtml='apptainer run $RIVET_CONTAINER_NAME rivet-mkhtml'
alias yodamerge='apptainer run $RIVET_CONTAINER_NAME yodamerge'
