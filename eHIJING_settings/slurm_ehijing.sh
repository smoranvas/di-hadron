END=100
for i in $(seq 1 $END); do
    for A in D C Fe Pb; do
	echo running slurm job for A=$A, i=$i
	sbatch ./run_ehijing.sh $A $i
    done
done
