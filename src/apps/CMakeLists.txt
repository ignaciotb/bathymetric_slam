project(slam_ws_apps)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

include_directories(include)

## APPS

add_executable(bathy_slam_real
       src/test_slam_real.cpp
)

target_link_libraries(bathy_slam_real
   gicp_registration
   utils_g2o
   submaps_tools
   graph_construction
   ceres_optimizer
   read_g2o
   data_transforms
   benchmark
   std_data
   slam_solver
)

add_executable(read_auv_data
       src/read_auv_data.cpp
)

target_link_libraries(read_auv_data
    std_data navi_data gsf_data all_data xtf_data csv_data data_transforms benchmark submaps_tools utils_visualization)

