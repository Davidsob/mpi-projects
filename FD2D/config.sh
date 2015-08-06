#!/bin/bash

rm -rf CMakeFiles *Cache.txt *.cmake

#MPI_PATH=/scr_davidson/davidson/contrib/mpich/bin
MPI_PATH=/usr/local/bin

cmake \
    -DCMAKE_CXX_COMPILER:PATH=$MPI_PATH/mpic++ \
    -DCMAKE_C_COMPILER:PATH=$MPI_PATH/mpicc \
    -DCMAKE_CXX_FLAGS:STRING='-Wall -fPIC -pipe -std=c++11 -stdlib=libc++' \
    -DCMAKE_C_FLAGS:STRING='-fPIC -pipe -Wall' \

make -j8
