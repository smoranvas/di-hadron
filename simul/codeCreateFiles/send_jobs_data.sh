#!/bin/bash
app=dataJob.slurm

start=1
option=$1
met_opt=${option:0:2}    

if [ "$met_opt" = "Ca" ] ; then 
	end=2630
elif [ "$met_opt" = "Fe" ] ; then 
	end=4000     
elif [ "$met_opt" = "Pb" ] ; then 
	end=3650
fi

for n in $(seq $start $end);do
    echo "sending $n";
    perc=$(($n*100/$end))
    echo "$perc % , for $met_opt"
    while [ -z "$(sbatch --export=opt=$option --array=$n-$n -J $met_opt-$n $app 2> /dev/null)" ];do
        echo "re-sending for $met_opt: $n";
    done;
done
printf "DONE!!\n"


