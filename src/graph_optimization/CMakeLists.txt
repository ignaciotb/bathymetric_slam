project(graph_optimization)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/g2o_tools)

include_directories(include)

# add_executable(create_sphere
#   src/create_sphere.cpp
# )

# Add libraries and executables
add_library(utils_g2o
    src/utils_g2o.cpp
)

add_library(graph_construction
    src/graph_construction.cpp
)

add_library(read_g2o
    src/read_g2o.cpp
)

add_library(ceres_optimizer
    src/ceres_optimizer.cpp
    src/types.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(utils_g2o PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

target_include_directories(graph_construction PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

target_include_directories(read_g2o PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

target_include_directories(ceres_optimizer PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

# set_target_properties(create_sphere PROPERTIES OUTPUT_NAME create_sphere${EXE_POSTFIX})

# target_link_libraries(create_sphere
# 	${G2O_LIBS}
# )

target_link_libraries(utils_g2o
	submaps_tools
	${G2O_LIBS}
)

target_link_libraries(graph_construction
	submaps_tools
	${G2O_LIBS}
)

target_link_libraries(ceres_optimizer
	ceres 
	${GFLAGS_LIBRARIES}
)

target_link_libraries(read_g2o 
	ceres 
	${GFLAGS_LIBRARIES}
)

# target_link_libraries(pose_graph_3d ceres ${GFLAGS_LIBRARIES})
