dir=/data/sebouh/eg2/python/
create_tuples(){
    #n=10000000
    n=10000000
    #n=1000000 #for testing

    #uncomment these lines to create the singles tuples
    
    #rm ${dir}/input_${1}.root
    #hadd ${dir}/input_${1}.root ${2}/*.pass2.root 
    #python3 makeSingleHadronTuple.py ${dir}/input_${1}.root:ntuple_data ${dir}/out_D_${1}.root:hadrons -N $n -T 1 &
    #python3 makeSingleHadronTuple.py ${dir}/input_${1}.root:ntuple_data ${dir}/out_${1}.root:hadrons -N $n -T 2 &
    #wait
    
    #args=--dX=.05 --dAngle=2
    #python3 ../../sidis_tuples/python/mix.py out_D_${1}.root mixed_D_${1}.root -N=$n -n=10 -e -s  --dX=.05 --dAngle=2 &
    #python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_${1}.root -N=$n -n=10 -e -s --dX=.05 --dAngle=2  &
    #python3 ~/sidis_tuples/python/mix.py out_D_${1}.root mixed_xQ2_D_${1}.root -N=$n -n=10 -e -s  --dX=.05 --dAngle=999 --dQ2=.3 &
    #python3 ~/sidis_tuples/python/mix.py out_${1}.root mixed_xQ2_${1}.root -N=$n -n=10 -e -s --dX=.05 --dAngle=999 --dQ2=.3 &

    python3 mix.py ${dir}/out_D_${1}.root ${dir}/mixed_no_ecuts_originalK_D_${1}.root -N=$n -n=10 -e -s  --dX=999 --dAngle=999 --dQ2=999 --original &
    python3 mix.py ${dir}out_${1}.root ${dir}/mixed_no_ecuts_originalK_${1}.root -N=$n -n=10 -e -s --dX=999 --dAngle=999 --dQ2=999 --original &
    #python3 ../../sidis_tuples/python/mix.py out_D_${1}.root mixed_no_ecuts_D_${1}.root -N=$n -n=10 -s &
    #python3 ../../sidis_tuples/python/mix.py out_${1}.root mixed_no_ecuts_${1}.root -N=$n -n=10 -s & 
    wait
}

create_tuples Fe /home/seba/di-hadron/data/fe/ &
create_tuples C /home/seba/di-hadron/data/ca/ &
create_tuples Pb /home/seba/di-hadron/data/pb/ &
wait

#for a in Fe Ca Pb; do
#    echo python compare_dphi.py mixed_${a}.root mixed_D_${a}.root
#    python3 compare_dphi.py mixed_${a}.root mixed_D_${a}.root
#done
#echo python3 compare_dphi.py mixed_Fe.root mixed_Pb.root
#python3 compare_dphi.py mixed_Fe.root mixed_Pb.root
#echo python3 compare_dphi.py mixed_Fe.root mixed_Ca.root
#python3 compare_dphi.py mixed_Fe.root mixed_Ca.root
