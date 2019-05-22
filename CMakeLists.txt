cmake_minimum_required(VERSION 3.0)

project(road_plane)


find_package(PCL 1.8 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)



INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})



add_executable(road_plane main.cpp)
target_link_libraries(road_plane ${PCL_LIBRARIES} Eigen3::Eigen)


