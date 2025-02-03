project(implant LANGUAGES CXX)
CPMFindPackage(
        NAME ViennaPS
        VERSION 3.1.0
        GIT_REPOSITORY "https://github.com/ViennaTools/ViennaPS")


# Use CPM.cmake to find GSL
CPMFindPackage(
    NAME GSL
    VERSION 2.7.1
    FIND_PACKAGE_ARGUMENTS "REQUIRED"
)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaCS)
target_link_libraries(${PROJECT_NAME} PUBLIC ViennaTools::ViennaPS)
target_link_libraries(${PROJECT_NAME} PUBLIC GSL::gsl GSL::gslcblas)

add_dependencies(ViennaCS_Examples ${PROJECT_NAME})
