cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ViennaCS_Python LANGUAGES CXX)

add_custom_target(${PROJECT_NAME} ALL)

# --------------------------------------------------------------------------------------------------------
# Global CMake Configuration
# └ As mentioned earlier we re-use the vtk python package libs, to do so we
#   set the rpath for our python modules to point to the vtkmodules folder.
# --------------------------------------------------------------------------------------------------------

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

if (NOT APPLE)
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN")
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../vtkmodules")
    list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../viennacs.libs")
else ()
    list(APPEND CMAKE_INSTALL_RPATH "@loader_path")
    list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../vtkmodules")
    list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../viennacs.libs")
endif ()

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("../cmake/cpm.cmake")

set(PYBIND11_FINDPYTHON ON)

CPMFindPackage(
        NAME pybind11
        VERSION 2.12.0
        GIT_REPOSITORY "https://github.com/pybind/pybind11")

# --------------------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------------------

set(VIENNACS_PYTHON_MODULE_NAME "viennacs")
set(VIENNACS_LIBRARY_OUTPUT_DIR ${CMAKE_BINARY_DIR})

if (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(VIENNACS_LIBRARY_OUTPUT_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    message(STATUS "[ViennaCS] Using Library Output Directory for bindings")
endif ()

# --------------------------------------------------------------------------------------------------------
# Setup Bindings
# --------------------------------------------------------------------------------------------------------

viennacore_setup_binding(viennacs2d ${VIENNACS_LIBRARY_OUTPUT_DIR})
viennacore_setup_binding(viennacs3d ${VIENNACS_LIBRARY_OUTPUT_DIR})

target_link_libraries(viennacs2d PRIVATE ViennaCS)
target_link_libraries(viennacs3d PRIVATE ViennaCS)

target_compile_definitions(viennacs2d PRIVATE -DVIENNACS_PYTHON_DIMENSION=2
        -DVIENNACS_MODULE_NAME=viennacs2d)
target_compile_definitions(viennacs3d PRIVATE -DVIENNACS_PYTHON_DIMENSION=3
        -DVIENNACS_MODULE_NAME=viennacs3d)

# --------------------------------------------------------------------------------------------------------
# Setup Lib-Folder
# --------------------------------------------------------------------------------------------------------

set(VIENNACS_LIB_FOLDER "${VIENNACS_LIBRARY_OUTPUT_DIR}/viennacs.libs")

# Not required for both targets, one will suffice
viennacore_setup_embree_env(viennacs2d "${VIENNACS_LIB_FOLDER}")
viennacore_setup_vtk_env(viennacs2d "${VIENNACS_LIB_FOLDER}")
viennacore_setup_tbb_env(viennacs2d "${VIENNACS_LIB_FOLDER}")

install(
        DIRECTORY "${VIENNACS_LIB_FOLDER}"
        DESTINATION .
        OPTIONAL)
