for a in C Fe Pb; do
    for i in 1 2 3 4 5 6 7 8 9; do
	python3 tmp2.py $a $i &
    done
done
wait
