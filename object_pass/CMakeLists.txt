project(Passer)
cmake_minimum_required(VERSION 2.8)

include_directories(
/scr_davidson/davidson/contrib/mpich/include/
)
# set the executable
add_executable(Passer main.cpp)
