project(implant LANGUAGES CXX)
CPMFindPackage(
        NAME ViennaPS
        VERSION 3.1.0
        GIT_REPOSITORY "https://github.com/ViennaTools/ViennaPS")


add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaCS)
target_link_libraries(${PROJECT_NAME} PUBLIC ViennaTools::ViennaPS)

add_dependencies(ViennaCS_Examples ${PROJECT_NAME})
