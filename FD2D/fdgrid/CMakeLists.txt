set(HEADERS
communicatorMap.h
LocalGrid.h
)

set(SOURCES
communicatorMap.cpp
LocalGrid.cpp
)

include_directories(
    ${PROJECT_SOURCE_DIR}/fdutils
    ${PROJECT_SOURCE_DIR}/boundaryConditions
)

set(EXTERNAL_LIBRARIES
fdutils
)

add_library(fdGrids ${HEADERS} ${SOURCES})
target_link_libraries(fdGrids ${EXTERNAL_LIBRARIES})
