project(libcerealization)

add_library(cerealization INTERFACE)


target_include_directories(cerealization INTERFACE
    #"${PROJECT_SOURCE_DIR}/include"
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

install(TARGETS cerealization EXPORT CerealizationConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.
install(EXPORT CerealizationConfig DESTINATION share/Cerealization/cmake)

# This makes the project importable from the build directory
export(TARGETS cerealization FILE CerealizationConfig.cmake)
