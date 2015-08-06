#!/bin/bash

rm -rf CMakeFiles *Cache.txt *.cmake

cmake \
    -DCMAKE_CXX_COMPILER:PATH=/scr_davidson/davidson/contrib/mpich/bin/mpic++ \
    -DCMAKE_C_COMPILER:PATH=/scr_davidson/davidson/contrib/mpich/bin//mpicc \
    -DCMAKE_CXX_FLAGS:STRING='-Wall -fPIC -pipe -std=c++11 -stdlib=libc++' \
    -DCMAKE_C_FLAGS:STRING='-fPIC -pipe -Wall' \
    -G Xcode
