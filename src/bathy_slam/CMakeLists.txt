project(bathy_slam)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(include)

# Add libraries and executables

add_library(slam_solver
    src/bathy_slam.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(slam_solver PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

target_link_libraries(slam_solver PUBLIC
    submaps_tools
    graph_construction
    gicp_registration
)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS slam_solver EXPORT BathySlamConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT BathySlamConfig DESTINATION share/BathySlam/cmake)

# This makes the project importable from the build directory
export(TARGETS slam_solver FILE BathySlamConfig.cmake)
