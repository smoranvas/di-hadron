create_tuples(){
    #n=10000000
    n=10000000

    #rm out_${1}.root
    #python3 makeSingleHadronTuple.py ${2}:ntuple_data out_${1}.root:hadrons -N $n -T 2 &
    wait
    #args=--dX=.05 --dAngle=2
    #python3 ../../sidis_tuples/python/mix.py out_D_${1}.root mixed_D_${1}.root -N=$n -n=10 -e -s  --dX=.05 --dAngle=2 &
    #python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_${1}.root -N=$n -n=10 -e -s --dX=.05 --dAngle=2  &
    
    python3 ../../sidis_tuples/python/mix.py ${2} mixed_xQ2_GiBUU_${1}.root -N=$n -n=10 -e -s --dX=.05 --dAngle=999 --dQ2=.3 &
    #python3 ../../sidis_tuples/python/mix.py out_D_${1}.root mixed_no_ecuts_D_${1}.root -N=$n -n=10 -s &
    #python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_no_ecuts_${1}.root -N=$n -n=10 -s & 
    wait
}

create_tuples D /home/sebouh/di-hadron/GiBUU_Singles_D.root &
create_tuples C /home/sebouh/di-hadron/GiBUU_Singles_C.root &
create_tuples Fe /home/sebouh/di-hadron/GiBUU_Singles_Fe.root &
create_tuples Pb /home/sebouh/di-hadron/GiBUU_Singles_Pb.root &
wait

#for a in Fe Ca Pb; do
#    echo python compare_dphi.py mixed_${a}.root mixed_D_${a}.root
#    python3 compare_dphi.py mixed_${a}.root mixed_D_${a}.root
#done
#echo python3 compare_dphi.py mixed_Fe.root mixed_Pb.root
#python3 compare_dphi.py mixed_Fe.root mixed_Pb.root
#echo python3 compare_dphi.py mixed_Fe.root mixed_Ca.root
#python3 compare_dphi.py mixed_Fe.root mixed_Ca.root
