cd ..; export PYTHIA8=`pwd`; cd -
rm -rf build/; mkdir build; cd build; cmake ..; make; cd ..
