create_tuples(){
    #rm input_${1}.root
    #hadd input_${1}.root ${2}/*.root 
    #python3 makeSingleHadronTuple.py input_${1}.root:ntuple_data out_D_${1}.root:hadrons -N 1000000 -T 1 &
    python3 makeSingleHadronTuple.py input_${1}.root:ntuple_sim out_${1}.root:hadrons -N 10000000 -T ${3} --isMC &
    wait
    #python3 ../../sidis_tuples/python/mix.py out_D_${1}.root mixed_D_${1}.root -N=10000000 -n=10 -e -s  --dX=.05 --dAngle=2 &
    #python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_${1}.root -N=100000000 -n=10 -e -s --dX=.05 --dAngle=2 --dQ2=9999  &
    python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_xQ2_${1}.root -N=100000000 -n=10 -e -s --dX=.05 --dAngle=999 --dQ2=.3  &
    #python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_no_ecuts_${1}.root -N=100000000 -n=10 -s --dX=999 --dAngle=999 --dQ2=999  &
    #wait
}

create_tuples MC_C /home/seba/di-hadron/simul/C/ 2 &
create_tuples MC_D2 /home/seba/di-hadron/simul/D2 1 &
create_tuples MC_Fe /home/seba/di-hadron/simul/Fe 2 &
create_tuples MC_Pb /home/seba/di-hadron/simul/Pb 2 &

wait


