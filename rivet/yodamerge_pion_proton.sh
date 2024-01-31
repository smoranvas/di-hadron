for A in D C Fe Pb; do
    cd /volatile/clas12/spaul/beagle_output/
    apptainer run $RIVET_CONTAINER_NAME yodamerge -o e${A}-beagle.yoda e${A}_*.yoda
done
