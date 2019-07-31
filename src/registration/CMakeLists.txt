project(registration)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(include)

# Add libraries and executables

# add_library(utils_registration
#     src/utils_registration.cpp
# )

add_library(utils_visualization
    src/utils_visualization.cpp
)

add_library(gicp_registration
    src/gicp_reg.cpp
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
# target_include_directories(utils_registration PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#     $<INSTALL_INTERFACE:include>
#     PRIVATE src)

target_include_directories(utils_visualization PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

target_include_directories(gicp_registration PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)


# Link the libraries
# target_link_libraries(utils_registration PUBLIC
#     submaps_tools    
#     ${PCL_LIBRARIES}
# )

target_link_libraries(utils_visualization PUBLIC
    submaps_tools
    graph_construction
    ceres_optimizer
    ${PCL_LIBRARIES}
)

target_link_libraries(gicp_registration PUBLIC
    utils_visualization
    submaps_tools
    benchmark
    ${PCL_LIBRARIES}
)
