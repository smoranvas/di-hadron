for a in D C Fe Pb ;do
    python3 GiBUU_PairNtuple.py $a tweak/GiBUU_${a}_tweak_harut.root GiBUU_Pairs_${a}_tweak_harut.root &
done
