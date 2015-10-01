project(FD2D)
cmake_minimum_required(VERSION 2.8)

add_subdirectory(fdgrid)
add_subdirectory(fdutils)
add_subdirectory(boundaryConditions)
add_subdirectory(fdmodels)
add_subdirectory(fddatamanager)

include_directories(
fdgrid
fdsources
fdutils
boundaryConditions
fdmodels
fddatamanager
/usr/local/include/
)

set(EXTERNAL_LIBRARIES
fdGrids
fdutils
boundaryConditions
models
fddata
)

# set the executable
#add_executable(Wave2D fd2d_wave.cpp)
#target_link_libraries(Wave2D ${EXTERNAL_LIBRARIES})

add_executable(Heat2D fd2d_heat.cpp)
target_link_libraries(Heat2D ${EXTERNAL_LIBRARIES})

add_executable(HeatAdvection2D fd2d_heat_advection.cpp)
target_link_libraries(HeatAdvection2D ${EXTERNAL_LIBRARIES})

add_executable(Transport2D fd2d_transport.cpp)
target_link_libraries(Transport2D ${EXTERNAL_LIBRARIES})

add_executable(Flow2D fd2d_compressible_flow.cpp)
target_link_libraries(Flow2D ${EXTERNAL_LIBRARIES})

add_executable(Shock2D fd2d_euler.cpp)
target_link_libraries(Shock2D ${EXTERNAL_LIBRARIES})
