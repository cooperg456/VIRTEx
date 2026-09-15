cmake_minimum_required(VERSION 3.25)

project(VIRTEx
    VERSION 0.0.1
    DESCRIPTION "Viral Infection, Recovery, and Transmission Explorer"
    LANGUAGES CXX CUDA
)

##  Nice try Microslop
if (WIN32)
    message(FATAL_ERROR "VIRTEx does not support Windows. Build using Linux, or WSL2 with GPU passthrough.")
endif ()

include(GNUInstallDirs)

################################################################################
##  general settings and options
################################################################################

##  VIRTEX_SHARED_LIBS option (undefined by default) can be used to force shared/static build
##  VIRTEX_CUDA_ARCH option can be set to native for faster compile times
option(VIRTEX_BUILD_EXECUTABLE "Build VIRTEx executable" ${PROJECT_IS_TOP_LEVEL})
option(VIRTEX_BUILD_TESTS "Build VIRTEx tests" OFF)
option(VIRTEX_BUILD_EXAMPLES "Build VIRTEx examples" OFF)
option(VIRTEX_BUILD_DOCS "Build VIRTEx documentation" OFF)
option(VIRTEX_INSTALL "Generate target for installing VIRTEx" ${PROJECT_IS_TOP_LEVEL})
set(VIRTEX_CUDA_ARCH all-major CACHE STRING "CUDA architectures")

if (DEFINED VIRTEX_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ${VIRTEX_SHARED_LIBS})
endif ()

if (PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

################################################################################
##  VIRTEx dependencies
################################################################################

find_package(CUDAToolkit REQUIRED)

################################################################################
##  VIRTEx library
################################################################################

add_library(virtex)
add_library(VIRTEx::virtex ALIAS virtex)

target_sources(virtex
    PRIVATE
        src/StochasticSimulation.cu
        src/FileIO.cpp
    PUBLIC
        FILE_SET HEADERS
            BASE_DIRS
                "${CMAKE_CURRENT_SOURCE_DIR}/include"
            FILES
                include/virtex/StochasticSimulation.hpp
                include/virtex/FileIO.hpp
)

target_compile_features(virtex
    PUBLIC cxx_std_20
    PRIVATE cuda_std_20)

target_link_libraries(virtex PRIVATE CUDA::curand)

set_target_properties(virtex PROPERTIES
    CXX_EXTENSIONS OFF
    CUDA_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CUDA_ARCHITECTURES "${VIRTEX_CUDA_ARCH}"
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_options(virtex PRIVATE
    "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:SHELL:-G -g>"
    "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Release,RelWithDebInfo>>:-O3>"
)

################################################################################
##  virtex executable
################################################################################

if (VIRTEX_BUILD_EXECUTABLE)
    add_executable(virtex-cli)
    add_executable(VIRTEx::virtex-cli ALIAS virtex-cli)

    target_sources(virtex-cli PRIVATE
        src/main.cpp
        src/CLITool.cpp
    )

    target_link_libraries(virtex-cli PRIVATE VIRTEx::virtex)

    set_target_properties(virtex-cli PROPERTIES
        OUTPUT_NAME virtex
        CXX_EXTENSIONS OFF
    )

    set_property(TARGET virtex-cli PROPERTY
        INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

################################################################################
##  install rules
################################################################################

if(VIRTEX_INSTALL AND NOT CMAKE_SKIP_INSTALL_RULES)
    install(TARGETS virtex
        LIBRARY
        ARCHIVE
        FILE_SET HEADERS
    )

    if(VIRTEX_BUILD_EXECUTABLE)
        install(TARGETS virtex-cli
            RUNTIME
        )
    endif()
endif()

################################################################################
##  other targets
################################################################################

if(VIRTEX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(VIRTEX_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(VIRTEX_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    doxygen_add_docs(docs include)
endif()
