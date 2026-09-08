cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(VIRTEx
    VERSION 0.0.1
    LANGUAGES CXX
)

option(INSTALL_MODELS "Install reaction network JSON files" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

##  virtex executable

add_executable(virtex)

target_sources(virtex PRIVATE
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/CLITool.cpp
    ${CMAKE_SOURCE_DIR}/src/FileIO.cpp
)

target_include_directories(virtex PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

##  platform dependant compilation

if(APPLE)

else()
    enable_language(CUDA)

    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES all-major)

    find_package(CUDAToolkit REQUIRED)

    set_target_properties(virtex PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )

    target_sources(virtex PRIVATE
        ${CMAKE_SOURCE_DIR}/src/StochasticSimulation.cu
    )

    target_link_libraries(virtex PRIVATE
        CUDA::curand
    )
endif()

##  install directories

include(GNUInstallDirs)

install(TARGETS virtex RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(INSTALL_MODELS)
    install(DIRECTORY models/
            DESTINATION ${CMAKE_INSTALL_DATADIR}/virtex/models
            FILES_MATCHING PATTERN "*.json")
endif()

##  build type

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug or Release)" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")

target_compile_options(virtex PRIVATE
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:-G -g>
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Release>>:-O3 --use_fast_math>
)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
