#!/bin/bash
#SBATCH --time=00:10:05
#SBATCH --mail-type=FAIL
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

met_opt=${opt:0:2}    ## this takes values: "Ca" or "Fe" or "Pb"


if [ "$met_opt" = "Ca" ] ; then
	rn_info="thickD2"
	met_lc="ca"  	
	Metal="C"
elif [ "$met_opt" = "Fe" ] ; then
	rn_info="thickD2"
	met_lc="fe" 
	Metal="Fe"
elif [ "$met_opt" = "Pb" ] ; then
	rn_info="thinD2"	
	met_lc="pb" 
	Metal="Pb"
fi

runlist="/user/s/smoran/jobs/listsAllRunNumbers/${Metal}.txt"

cd $TMPDIR
exe="GetSimpleTuple_data"
sourceDir="/user/s/smoran/di-hadron_ANA/GetSimpleTuple/bin"
sourceDir2="/user/s/smoran/di-hadron_ANA/GetSimpleTuple/src"
outdir="/work/smoran/di-hadron/data/P25GeV"

cp $sourceDir/$exe exe
if [ ! -d $outdir/$met_lc ]; then mkdir -p $outdir/$met_lc; fi
    
rn=`sed -n "$SLURM_ARRAY_TASK_ID{p;q}" $runlist` 

find /data/jlab/mss/clas/eg2a/production/Pass2/Clas/ -name $rn -exec cp {} in.root \; 
find . -name "*.root" -print > dataFiles.txt
cat dataFiles.txt

./exe -t$Metal -r$rn
cp  "out.root" $outdir/$met_lc"/N_"$rn
cp  /user/s/smoran/di-hadron_ANA/GetSimpleTuple/include/TIdentificatorV2.h $outdir/$met_lc"/code_used.h"

