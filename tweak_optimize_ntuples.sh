for b in `ls tweak_optimize/2023_02_10*.root`; do
    b=`basename $b`
    for a in D  ;do
	python3 GiBUU_PairNtuple.py $a tweak_optimize/$b GiBUU_Pairs_tweak_${b} &
    done
done
