cmake_minimum_required (VERSION 2.6)
project (MeshFix)
include_directories(JMeshExt-1.0alpha_src/include ../JMeshLib-1.1/include JMeshExt-1.0alpha_src/OpenNL3.0.2/src/NL JMeshExt-1.0alpha_src/OpenNL3.0.2/src)
link_directories(${LINK_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}/JMeshExt-1.0alpha_src/lib)
link_directories(${LINK_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}/JMeshExt-1.0alpha_src/OpenNL3.0.2/binaries/lib)
set(meshfix_h 
JMeshExt-1.0alpha_src/include/jrs_predicates.h
JMeshExt-1.0alpha_src/include/detectIntersections.h
JMeshExt-1.0alpha_src/include/exttrimesh.h
JMeshExt-1.0alpha_src/include/holeFilling.h
JMeshExt-1.0alpha_src/include/sparseLSystem.h
)
set(meshfix_src
JMeshExt-1.0alpha_src/src/jrs_predicates.c
JMeshExt-1.0alpha_src/src/detectIntersections.cpp
JMeshExt-1.0alpha_src/src/holeFilling.cpp
JMeshExt-1.0alpha_src/src/sparseLSystem.cpp
)
add_executable(meshfix ${meshfix_src} ${meshfix_h} meshfix.cpp)
target_link_libraries(meshfix jmesh nl)

