set(HEADERS
boundaryCondition.h
dirchletBoundaryCondition.h
neumannBoundaryCondition.h
convectiveCooling.h
)

set(SOURCES
dirchletBoundaryCondition.cpp
neumannBoundaryCondition.cpp
convectiveCooling.cpp
)

add_library(boundaryConditions ${HEADERS} ${SOURCES})
