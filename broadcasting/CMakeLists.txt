project(mpi_ring)
cmake_minimum_required(VERSION 2.8)

# set the executable
add_executable(BCast broadcast.cpp)
add_executable(Compare bcast_compare.cpp)
