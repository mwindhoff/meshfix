cmake_minimum_required (VERSION 2.6)
project (MeshFix)
add_subdirectory(contrib/JMeshLib)
include_directories(include contrib/JMeshLib/include contrib/OpenNL3.0.2/src/NL contrib/OpenNL3.0.2/src contrib/octree)
link_directories(${LINK_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}/contrib/JMeshLib/lib)
link_directories(${LINK_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}/contrib/OpenNL3.0.2/binaries/lib)
set(meshfix_h 
include/jrs_predicates.h
include/detectIntersections.h
include/exttrimesh.h
include/holeFilling.h
include/sparseLSystem.h
include/triangleoctree.h
)
set(meshfix_src
src/jrs_predicates.c
src/detectIntersections.cpp
src/holeFilling.cpp
src/sparseLSystem.cpp
src/mwExtensions.cpp
src/triangleoctree.cpp
)
add_custom_command(
	TARGET meshfix
	PRE_LINK
	COMMAND make
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
contrib/JMeshLib/
)
add_executable(meshfix ${meshfix_src} ${meshfix_h} meshfix.cpp)
target_link_libraries(meshfix jmesh nl)
