for a in D C Fe Pb ;do
    python3 GiBUU_PairNtuple.py $a tweak/${a}_tweak.root GiBUU_Pairs_${a}_tweak.root &
done
