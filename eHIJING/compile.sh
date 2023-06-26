rm -rf build
mkdir build
cd build
cmake .. -DPYTHIA8_LIBRARY=~/eHIJING-pythia/lib/libpythia8.a -DPYTHIA8_INCLUDE_DIR=~/eHIJING-pythia/include
make
cd ..
