create_tuples(){
    #rm input_${1}.root
    #hadd input_${1}.root ${2}/*.pass2.root 
    #python3 makeSingleHadronTuple.py input_${1}.root:ntuple_data out_D_${1}.root:hadrons -N 1000000 -T 1 &
    #python3 makeSingleHadronTuple.py input_${1}.root:ntuple_data out_${1}.root:hadrons -N 1000000 -T 2 &
    #wait
    #args=--dX=.05 --dAngle=2
    python3 ../../sidis_tuples/python/mix.py out_D_${1}.root mixed_no_ecuts_D_${1}.root -N=10000000 -n=10 -s  --dX=.05 --dAngle=2 &
    python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_no_ecuts_${1}.root -N=10000000 -n=10 -s --dX=.05 --dAngle=2  &
    wait
}

create_tuples Fe /home/seba/di-hadron/data/fe/ &
create_tuples C /home/seba/di-hadron/data/ca/ &
create_tuples Pb /home/seba/di-hadron/data/pb/ &
wait

for a in Fe Ca Pb; do
    echo python compare_dphi.py mixed_${a}.root mixed_D_${a}.root
    python3 compare_dphi.py mixed_${a}.root mixed_D_${a}.root
done
echo python3 compare_dphi.py mixed_Fe.root mixed_Pb.root
python3 compare_dphi.py mixed_Fe.root mixed_Pb.root
echo python3 compare_dphi.py mixed_Fe.root mixed_Ca.root
python3 compare_dphi.py mixed_Fe.root mixed_Ca.root
