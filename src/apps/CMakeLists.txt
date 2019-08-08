project(slam_ws_apps)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

include_directories(include)

## TESTS
add_executable(test_slam_simulation
        src/test_graph_optimization.cpp
)

add_executable(test_slam_real
       src/test_slam_real.cpp
)

target_link_libraries(test_slam_simulation
    gicp_registration
    utils_g2o
    submaps_tools
    graph_construction
    ceres_optimizer
    read_g2o
    data_transforms
    benchmark
)

target_link_libraries(test_slam_real
   gicp_registration
   utils_g2o
   submaps_tools
   graph_construction
   ceres_optimizer
   read_g2o
   benchmark
   std_data
)

