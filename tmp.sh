for a in Fe Pb ;do
    python3 GiBUU_PairNtuple.py $a tweak/GiBUU_${a}_tweak_harut_fullstats.root GiBUU_Pairs_${a}_tweak_harut.root &
done
