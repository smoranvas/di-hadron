for t in C Fe Pb D2_pb;
do 
    mkdir /work/sebouh/di-hadron/simul/${t};
    for n in 1 2 3 4 5 6 7 8 9; 
    do 
	hadd_out=/work/sebouh/di-hadron/simul/${t}/${n}.root
	rm $hadd_out
	hadd -k $hadd_out /work/sebouh/di-hadron/simul/${t}${n}*/ROOT/*.root; 
    done 
    

done
