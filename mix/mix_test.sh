dir=/data/sebouh/eg2/python/

#edit the tag and args for each variation of the event-mixing
#tag=no_ecuts_originalK
#args="--dX=999 --dAngle=999 --dQ2=999 --original"

tag=angle_cut_5deg
args="-e --dAngle=5 --dX=999 --dQ2=9999" 

tag=no_ecuts_originalPT
args="--dX=999 --dAngle=999 --dQ2=999 --originalPT"  

echo args are $args

create_tuples(){

    n=1000000 #for testing
    #uncomment these lines to create the singles tuples
    
    #rm ${dir}/input_${1}.root
    #hadd ${dir}/input_${1}.root ${2}/*.pass2.root 
    #python3 makeSingleHadronTuple.py ${dir}/input_${1}.root:ntuple_data ${dir}/out_D_${1}.root:hadrons -N $n -T 1 &
    #python3 makeSingleHadronTuple.py ${dir}/input_${1}.root:ntuple_data ${dir}/out_${1}.root:hadrons -N $n -T 2 &
    #wait

    python3 mix.py ${dir}/out_${1}.root ${dir}/mixed_${tag}_${1}.root -N=$n -n=10 -s  $args &
    
    python3 mix.py ${dir}/out_D_${1}.root ${dir}/mixed_${tag}_D_${1}.root -N=$n -n=10 -s  $args &
    
}



create_tuples Fe /home/seba/di-hadron/data/fe/ 
create_tuples C /home/seba/di-hadron/data/ca/  
create_tuples Pb /home/seba/di-hadron/data/pb/ 

wait

hadd -f ${dir}/mixed_${tag}_D.root ${dir}/mixed_${tag}_D_*.root
