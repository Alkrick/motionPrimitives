cmake_minimum_required(VERSION 3.0.2)
project(motionPrimitives)

find_package(Eigen3 REQUIRED)

list(APPEND MP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND MP_INCLUDE_DIRS /include/)
set(MP_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(ProMP)