project(randomwalk)
cmake_minimum_required(VERSION 2.8)

add_subdirectory(randomnumbers)
include_directories(
    randomnumbers
    /scr_davidson/davidson/contrib/mpich/include/
)
set(EXTERNAL_LIBS
    Random
)
# set the executable
add_executable(RandomWalk random_walk.cpp)
target_link_libraries(RandomWalk ${EXTERNAL_LIBS})
