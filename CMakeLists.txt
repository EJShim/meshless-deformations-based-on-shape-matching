cmake_minimum_required(VERSION 2.8)

PROJECT(physics)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)



set(SRC_FILES
main
Mainwindow
Beam
PickInteractor
)

add_executable(physics ${SRC_FILES})


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
target_link_libraries(physics ${VTK_LIBRARIES})