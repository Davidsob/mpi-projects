#!/bin/bash

rm -rf CMakeFiles *Cache.txt *.cmake

#MPI_PATH=/scr_davidson/davidson/contrib/mpich/bin
MPI_PATH=/opt/local/bin

cmake \
    -DCMAKE_CXX_COMPILER:PATH=$MPI_PATH/mpicxx \
    -DCMAKE_C_COMPILER:PATH=$MPI_PATH/mpicc \
    -DCMAKE_CXX_FLAGS:STRING='-config=mpiconfig -Wall -fPIC -pipe -std=c++11 -Wno-overloaded-virtual' \
    -DCMAKE_C_FLAGS:STRING='-cc=gcc -fPIC -pipe -Wall -Wno-overloaded-virtual' \

make -j8
