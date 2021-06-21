#!/bin/bash
app=simJob.slurm

start=1
option=$1
if [ "$1" = "C" ] ; then 
	end=4439  
	met_opt="C"
elif [ "$1" = "Fe" ] ; then 
	end=4200  
	met_opt="Fe"
elif [ "$1" = "Pb" ]; then 
	end=4434  
	met_opt="Pb"
elif [ "$1" = "D2" ]; then 
	end=4440  
	met_opt="D2"
fi

for n in $(seq $start $end);do
    echo "sending $n";
    perc=$(($n*100/$end))
    echo "$perc % , for $met_opt"
    while [ -z "$(sbatch --export=opt=$option --array=$n-$n -J $met_opt-$n $app 2> /dev/null)" ];do
        echo "re-sending $n";
    done;
done
printf "DONE!!\n"

