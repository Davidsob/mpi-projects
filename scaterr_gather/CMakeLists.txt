project(mpi_ring)
cmake_minimum_required(VERSION 2.8)

add_subdirectory(randomnumbers)
include_directories(
    randomnumbers
)
# set the executable
add_executable(Avg mpi_avg.cpp)
target_link_libraries(Avg Random)

