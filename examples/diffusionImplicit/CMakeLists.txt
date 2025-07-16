find_package (Eigen3 3.3 NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
  message(STATUS "[ViennaCS] Eigen library not found, skipping example")
  return()
endif()


project(diffusionImplicit LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaCS)
configure_file(config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)

add_dependencies(ViennaCS_Examples ${PROJECT_NAME})
