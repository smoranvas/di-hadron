#!/bin/bash
#SBATCH --time=00:07:05
#SBATCH --mail-user=sebastian.moran@alumnos.usm.cl
#SBATCH --mail-type=FAIL
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

cd $TMPDIR
echo working on `pwd`

if [ "$opt" = "C" ] ; then
	list="sim_hayk_rn_ca.txt"
	met="C"
elif [ "$opt" = "Fe" ] ; then
	list="sim_hayk_rn_fe.txt"
	met="Fe"
elif [ "$opt" = "Pb" ]; then
	list="sim_hayk_rn_pb.txt"
	met="Pb"
elif [ "$opt" = "D2" ]; then
	list="sim_hayk_rn_d2.txt"
	met="D"
fi
cp /user/s/smoran/jobs/simul_jobs/old_list_format/$list list.txt
dfin=`sed -n "$SLURM_ARRAY_TASK_ID{p;q}" list.txt`
dfin=$dfin
dir=`dirname $dfin`
fin=`basename $dfin`
printf "copying the source file\n"
cp /eos/user/o/orsosa/HSim/$dfin in.root
printf "Creating the final directories \n"
outdir="/work/smoran/di-hadron/simul/$dir"
printf "Creating the executables variables\n"

exe="GetSimpleTuple_sim"

printf "Copying the executable \n"
sourceDir="/user/s/smoran/di-hadron_ANA/GetSimpleTuple/bin"
sourceDir2="/user/s/smoran/di-hadron_ANA/GetSimpleTuple/src"

cp $sourceDir/$exe exe
if [ ! -d $outdir ]; then mkdir -p $outdir; fi
fout=$outdir/$fin
printf "Executing...\n"
./exe -t$met -r1
printf "Copying the results to their final destination \n"
cp out.root $fout
printf "Copying the codes used in each case in the final destination \n"
cp $sourceDir2/$exe.cxx $outdir/code_used.cxx
printf "Done!"
