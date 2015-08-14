set(HEADERS
FDModel.h
FDHeatTransfer.h
FDTransport.h
FDCompressibleFlow.h
)

set(SOURCES
FDModel.cpp
FDHeatTransfer.cpp
FDTransport.cpp
FDCompressibleFlow.cpp
)

include_directories(
${PROJECT_SOURCE_DIR}/fdutils
${PROJECT_SOURCE_DIR}/boundaryConditions
${PROJECT_SOURCE_DIR}/fdgrid
${PROJECT_SOURCE_DIR}/fdsources
)

set(EXTERNAL_LIBS
    fdGrids
    fdutils
    boundaryConditions
)

add_library(models ${HEADERS} ${SOURCES})
target_link_libraries(models ${EXTERNAL_LIBS})
