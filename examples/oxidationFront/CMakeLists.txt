find_package(Eigen3 3.3.2 NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
  message(STATUS "[ViennaCS] Eigen library not found, skipping example")
  return()
endif()


project(oxidationFront LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaCS Eigen3::Eigen)
configure_file(config.txt ${CMAKE_CURRENT_BINARY_DIR}/config.txt COPYONLY)
configure_file(Efield.csv ${CMAKE_CURRENT_BINARY_DIR}/Efield.csv COPYONLY)
configure_file(writeCSV.py ${CMAKE_CURRENT_BINARY_DIR}/writeCSV.py COPYONLY)

add_dependencies(ViennaCS_Examples ${PROJECT_NAME})
