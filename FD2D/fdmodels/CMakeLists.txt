set(HEADERS
FDModel.h
FDHeatTransfer.h
FDTransport.h
FDCompressibleFlow.h
FDCompressibleEuler.h
)

set(SOURCES
FDModel.cpp
FDHeatTransfer.cpp
FDTransport.cpp
FDCompressibleFlow.cpp
FDCompressibleEuler.cpp
)

include_directories(
${PROJECT_SOURCE_DIR}/fdutils
${PROJECT_SOURCE_DIR}/boundaryConditions
${PROJECT_SOURCE_DIR}/fdgrid
${PROJECT_SOURCE_DIR}/fdsources
${PROJECT_SOURCE_DIR}/fddatamanager
)

set(EXTERNAL_LIBS
    fdGrids
    fdutils
    boundaryConditions
    fddata
)

add_library(models ${HEADERS} ${SOURCES})
target_link_libraries(models ${EXTERNAL_LIBS})
