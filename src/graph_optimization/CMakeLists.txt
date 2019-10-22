project(graph_optimization)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/g2o_tools)

#include_directories(include)

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

target_link_libraries(utils_g2o
	submaps_tools
	${G2O_LIBS}
)

target_link_libraries(graph_construction
	submaps_tools
        utils_g2o
	${G2O_LIBS}
)

target_link_libraries(ceres_optimizer
	ceres 
        submaps_tools
	${GFLAGS_LIBRARIES}
)

target_link_libraries(read_g2o 
	ceres 
	${GFLAGS_LIBRARIES}
)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS utils_g2o graph_construction ceres_optimizer read_g2o EXPORT GraphOptimizationConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT GraphOptimizationConfig DESTINATION share/GraphOptimization/cmake)

# This makes the project importable from the build directory
export(TARGETS utils_g2o graph_construction ceres_optimizer read_g2o FILE GraphOptimizationConfig.cmake)
