#!/bin/bash
echo "running job"
cd /home/spaul/di-hadron/beagle
source ./SDCC_setup.sh
source ./beagle_setup.sh
echo "sourced setup"
A=$1
i=$2
echo A=${A} i=${i}

output_dir=/volatile/clas12/spaul/beagle_output
log_dir=/volatile/clas12/spaul/logs/
mkdir -p $output_dir
mkdir -p $log_dir

txt_file=${output_dir}/e${A}_${i}.txt

#copy the input files and make versions of them with the correct output file names, and use different seeds in Pythia
#first create a temp folder (t) for to put the config files in
mkdir -p t
cat eAS1noq | sed 's|OUTFILE.txt|'$txt_file'|g' | sed 's|1234567|12345'${i}'|g' > t/${A}_${i}
cat input/e${A}_clas.inp | sed 's|e'${A}'S1noq|t/'${A}'_'${i}'|g' > t/e${A}_${i}_clas.inp
$BEAGLESYS/BeAGLE < t/e${A}_${i}_clas.inp > ${log_dir}/e${A}_${i}_clas.log

#convert file to root, then hepmc3
root_file=${output_dir}'/'`basename -s .txt ${txt_file}`'.root'
echo 'BuildTree("'${txt_file}'","'${output_dir}'")' | eic-smear
echo 'TreeToHepMC("'${root_file}'","'${output_dir}'")' | eic-smear

##comment out these lines if you don't want to save the large text file and the root file
rm $root_file
rm $txt_file
