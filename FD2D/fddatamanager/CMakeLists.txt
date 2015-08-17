set(HEADERS
DataManager.h
FDDataManager.h
)

set(SOURCES
FDDataManager.cpp
)

include_directories(
${PROJECT_SOURCE_DIR}/fdutils
${PROJECT_SOURCE_DIR}/fdsources
)

set(EXTERNAL_LIBS
    fdutils
)

add_library(fddata ${HEADERS} ${SOURCES})
target_link_libraries(fddata ${EXTERNAL_LIBS})
