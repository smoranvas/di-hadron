#!/bin/bash
#export PYTHIA8DATA=/global/homes/w/wk42/miniconda3/envs/ehijing/share/Pythia8/xmldoc
#export PYTHIA8=/global/homes/w/wk42/miniconda3/envs/ehijing
label=$1
i=$2

if [[ $label == "D" ]]; then
    A=2
    Z=1
    echo "deuterium"
elif [[ $label == "C" ]]; then
    A=12
    Z=6
    echo "carbon"
elif [[ $label == "Fe" ]]; then
    A=56
    Z=26
    echo "iron"
elif [[ $label == "Pb" ]]; then
    A=207
    Z=82
    echo "lead"
else 
    echo "first argument must be D, C, Fe, or Pb"
    exit
fi


exe=build/ehijing-test
#Neve=10000000
Neve=1000000
#Configfile=s20.setting
#Configfile=harut.setting
#Configfile2=hadronization.setting
#topLevelFolder=EventsTweak

Configfile=ahmed.setting
Configfile2=hadronization_ahmed.setting
topLevelFolder=/volatile/clas12/spaul/ehijing_output/



K=4.0
M=1 # Generlizaed HT:1,  HIgher-twist:0, both in the soft gluon emission limit.


folder=$topLevelFolder/$label
TablePath="Tables/$K"
mkdir -p $folder
mkdir -p $TablePath

$exe $Neve $Z $A $M $K $TablePath $folder $Configfile $Configfile2 > $folder/log_${i}.txt 
