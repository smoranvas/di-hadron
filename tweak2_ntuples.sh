for a in D C Fe Pb ;do
    python3 GiBUU_PairNtuple.py $a tweak/${a}_tweak2.root GiBUU_Pairs_${a}_tweak2.root &
done
