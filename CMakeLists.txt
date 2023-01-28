cmake_minimum_required(VERSION 3.0.2)
project(motionPrimitives)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIR}
)

add_executable(ProMP_test src/ProMP/ProMP_Test.cpp
src/ProMP/ProMP.cpp
)

target_link_libraries(ProMP_test 
Eigen3::Eigen
)

add_library(promp src/ProMP/ProMP.cpp)

install(TARGETS promp DESTINATION /usr/lib)