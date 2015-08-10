set(HEADERS
boundaryCondition.h
dirchletBoundaryCondition.h
neumannBoundaryCondition.h
convectiveCooling.h
boundaryConditions.h
)

set(SOURCES
dirchletBoundaryCondition.cpp
neumannBoundaryCondition.cpp
convectiveCooling.cpp
)

add_library(boundaryConditions ${HEADERS} ${SOURCES})
