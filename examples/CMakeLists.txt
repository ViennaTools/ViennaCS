add_custom_target(ViennaCS_Examples ALL)

set(VIENNACS_ARTIFACTS_DIRECTORY $<1:${PROJECT_BINARY_DIR}/artifacts>)

if(WIN32)
  viennacore_setup_embree_env(ViennaCS_Examples ${VIENNACS_ARTIFACTS_DIRECTORY})
  viennacore_setup_vtk_env(ViennaCS_Examples ${VIENNACS_ARTIFACTS_DIRECTORY})
  viennacore_setup_tbb_env(ViennaCS_Examples ${VIENNACS_ARTIFACTS_DIRECTORY})
endif()

viennacore_add_subdirs(${CMAKE_CURRENT_LIST_DIR})
