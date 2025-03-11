cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(
        ViennaCS
        LANGUAGES CXX
        VERSION 1.0.0)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

option(VIENNACS_BUILD_PYTHON "Build python bindings" ON)

option(VIENNACS_BUILD_EXAMPLES "Build examples" ON)
option(VIENNACS_BUILD_TESTS "Build tests" ON)

set(VIENNACS_LOOKUP_DIRS
        ""
        CACHE STRING "Directories to account for when searching installed dependencies")

list(APPEND CMAKE_PREFIX_PATH ${VIENNACS_LOOKUP_DIRS})

# --------------------------------------------------------------------------------------------------------
# Setup Library
# --------------------------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} INTERFACE)
add_library(ViennaTools::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES CXX_STANDARD 17
        CXX_EXTENSIONS OFF
        CXX_STANDARD_REQUIRED ON)

if (MSVC)
    # https://learn.microsoft.com/cpp/c-runtime-library/math-constants
    # TODO: In case C++20 is adpoted any time soon: https://cppreference.com/w/cpp/numeric/constants
    target_compile_definitions(${PROJECT_NAME} INTERFACE _USE_MATH_DEFINES=1)
endif ()

# --------------------------------------------------------------------------------------------------------
# Include directories
# --------------------------------------------------------------------------------------------------------

target_include_directories(
        ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/viennacs>
        $<INSTALL_INTERFACE:include/viennacs-${PROJECT_VERSION}>)

# --------------------------------------------------------------------------------------------------------
# Check whether or not dependencies are installed on the system
# └ VTK requires C Language support (HDF5)
# --------------------------------------------------------------------------------------------------------

enable_language(C)

find_package(ViennaRay QUIET)
find_package(ViennaLS QUIET)

set(VIENNACS_SYSTEM_VIENNARAY
        ${ViennaRay_FOUND}
        CACHE INTERNAL "")

set(VIENNACS_SYSTEM_VIENNALS
        ${ViennaLS_FOUND}
        CACHE INTERNAL "")

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("cmake/cpm.cmake")

CPMAddPackage(
        NAME ViennaCore
        VERSION 1.3.0
        GIT_REPOSITORY "https://github.com/ViennaTools/ViennaCore"
        OPTIONS "VIENNACORE_FORMAT_EXCLUDE docs/")

CPMAddPackage(
        NAME PackageProject
        VERSION 1.13.0
        GIT_REPOSITORY "https://github.com/TheLartians/PackageProject.cmake"
        EXCLUDE_FROM_ALL ${VIENNACS_BUILD_PYTHON})

CPMFindPackage(
        NAME ViennaRay
        VERSION 3.3.0
        GIT_REPOSITORY "https://github.com/ViennaTools/ViennaRay"
        EXCLUDE_FROM_ALL ${VIENNACS_BUILD_PYTHON})

CPMFindPackage(
        NAME ViennaLS
        VERSION 4.2.2
        GIT_REPOSITORY "https://github.com/ViennaTools/ViennaLS"
        EXCLUDE_FROM_ALL ${VIENNACS_BUILD_PYTHON})

target_link_libraries(${PROJECT_NAME} INTERFACE ViennaTools::ViennaCore ViennaTools::ViennaLS
        ViennaTools::ViennaRay)

# --------------------------------------------------------------------------------------------------------
# Setup Examples
# --------------------------------------------------------------------------------------------------------

if (VIENNACS_BUILD_EXAMPLES)
    message(STATUS "[ViennaCS] Building Examples")
    add_subdirectory(examples)
endif ()

# --------------------------------------------------------------------------------------------------------
# Setup Tests
# --------------------------------------------------------------------------------------------------------

if (VIENNACS_BUILD_TESTS)
    message(STATUS "[ViennaCS] Building Tests")

    enable_testing()
    add_subdirectory(tests)
endif ()

# --------------------------------------------------------------------------------------------------------
# Setup Python Bindings
# --------------------------------------------------------------------------------------------------------

if (VIENNACS_BUILD_PYTHON)
    message(STATUS "[ViennaPS] Building Python Bindings")
    add_subdirectory(python)
    return()
endif ()

# --------------------------------------------------------------------------------------------------------
# Install Target
# --------------------------------------------------------------------------------------------------------

packageProject(
        NAME ${PROJECT_NAME} NAMESPACE ViennaTools
        VERSION ${PROJECT_VERSION}
        BINARY_DIR ${PROJECT_BINARY_DIR}
        INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/viennacs
        INCLUDE_DESTINATION include/viennacs-${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
        DEPENDENCIES "ViennaCore;ViennaLS;ViennaRay")
