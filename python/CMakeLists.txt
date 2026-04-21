cmake_minimum_required(VERSION 3.14)
project(ViennaCS_Python LANGUAGES CXX)

add_custom_target(${PROJECT_NAME} ALL)

# --------------------------------------------------------------------------------------------------------
# Global CMake Configuration
# └ See: https://github.com/ViennaTools/ViennaPS/blob/c76e371817a797dfe2800691f00cb93317b731fa/python/CMakeLists.txt#L8
# --------------------------------------------------------------------------------------------------------

set(CMAKE_MACOSX_RPATH ON)

set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

if(NOT APPLE)
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN")
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../vtkmodules")
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../viennacs.libs")
else()
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path")
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../vtkmodules")
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../viennacs.libs")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("../cmake/cpm.cmake")

set(PYBIND11_FINDPYTHON ON)

CPMFindPackage(
  NAME pybind11
  VERSION 3.0.0
  GIT_REPOSITORY "https://github.com/pybind/pybind11")

# --------------------------------------------------------------------------------------------------------
# Setup Bindings
# --------------------------------------------------------------------------------------------------------

set(PKG viennacs)
set(MOD _core)

pybind11_add_module(${MOD} pyWrap.cpp)
target_link_libraries(${MOD} PRIVATE ViennaCS)
target_compile_definitions(${MOD} PRIVATE VIENNACS_MODULE_NAME=${MOD})

# Put the built .so in build/<pkg> for convenience
set_target_properties(
  ${MOD}
  PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
             RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
             ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}")

# Install into site-packages/<pkg>/  (RELATIVE path!)
install(
  TARGETS ${MOD}
  LIBRARY DESTINATION ${PKG}
  RUNTIME DESTINATION ${PKG}
  ARCHIVE DESTINATION ${PKG})

# Package files
install(FILES ${CMAKE_SOURCE_DIR}/python/__init__.py DESTINATION ${PKG})
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/data/implant
  DESTINATION ${PKG}/data
  OPTIONAL)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/data/ImplData/Taurus
  DESTINATION ${PKG}/data/ImplData
  OPTIONAL)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/data/AnnealData
  DESTINATION ${PKG}/data
  OPTIONAL)

add_dependencies(${PROJECT_NAME} ${MOD})

# --------------------------------------------------------------------------------------------------------
# Setup Lib-Folder
# --------------------------------------------------------------------------------------------------------

set(VIENNACS_LIB_FOLDER ${CMAKE_BINARY_DIR}/${PKG}.libs)

viennacore_setup_vtk_env(${MOD} ${VIENNACS_LIB_FOLDER})

install(
  DIRECTORY ${VIENNACS_LIB_FOLDER}
  DESTINATION .
  OPTIONAL)
