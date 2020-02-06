project(meas_models)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/submaps)

include_directories(include)

add_library(multibeam_simple
    src/multibeam_simple.cpp
)

#Define headers for this library. PUBLIC headers are used for
#compiling the library, and will be added to consumers' build
#paths.
target_include_directories(multibeam_simple PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

target_link_libraries(multibeam_simple PUBLIC
    submaps_tools
    ${PCL_LIBRARIES}
)
