cmake_minimum_required(VERSION 3.0)

project(road_plane)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

add_executable(road_plane main.cpp)
target_link_libraries(road_plane)


