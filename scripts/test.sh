#!/bin/bash
rm -rf output/hdf5*
make 
for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
echo ./src/falle input/ryujones${i} 
./src/falle input/ryujones${i} > /dev/null
cd output 
#ls -t hdf5* | head -1
cp $(ls -t out_2d_* | head -1) ryujones${i}.txt
#name=test${i}.h5
#mv `ls -t hdf5* | head -1` name
cd ..
done
