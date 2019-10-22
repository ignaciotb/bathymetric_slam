cmake_minimum_required(VERSION 2.8)
project(uw_slam_projects)

# check c++14 / c++0x
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -fPIC")
     add_compile_options(-std=c++14 -Wall -Wextra)
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fPIC")
     add_compile_options(-std=c++11 -Wall -Wextra)
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 or C++11 support. Please use a different C++ compiler.")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -fPIC")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set( CMAKE_BUILD_TYPE "Release" )

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL 1.8 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(g2o REQUIRED)

# AUVLIB stuff
 # find_package(LIBIGL REQUIRED QUIET)
 find_package(LibXTF REQUIRED)
 find_package(LibGSF REQUIRED)
 find_package(EigenCereal REQUIRED)
 find_package(DataTools REQUIRED)
 find_package(SonarTracing REQUIRED)
# find_package(AuvlibGlad REQUIRED)
 find_package(BathyMaps REQUIRED)

link_directories(${PCL_LIBRARY_DIRS})
link_directories(${G2O_LIBRARY_DIRS})

add_definitions(${PCL_DEFINITIONS})
add_definitions(-DCERES_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})

include_directories(
	${CERES_INCLUDE_DIRS}
	${EIGEN_INCLUDE_DIRS} 
	${G2O_INCLUDE_DIRS} 
	${PCL_INCLUDE_DIRS}
)

LIST(APPEND G2O_LIBS
  cxsparse
  cholmod
  g2o_cli g2o_ext_freeglut_minimal g2o_simulator
  g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d
  g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only
  g2o_types_sba g2o_types_slam3d g2o_csparse_extension
  g2o_opengl_helper g2o_solver_dense g2o_stuff
  g2o_types_sclam2d g2o_parser g2o_solver_pcg
  g2o_types_data g2o_types_sim3
)

add_subdirectory(apps)
add_subdirectory(graph_optimization)
add_subdirectory(registration)
add_subdirectory(submaps_tools)
add_subdirectory(bathy_slam)
