project(mpi_ring)
cmake_minimum_required(VERSION 2.8)

add_subdirectory(randomnumbers)
include_directories(
    randomnumbers
)
set(EXTERNAL_LIBS
    Random
)
# set the executable
add_executable(CheckStatus check_status.cpp)
target_link_libraries(CheckStatus ${EXTERNAL_LIBS})
