add_executable(ProMP_LearnEx ProMP_LearnEx.cpp
ProMP.cpp
)

target_include_directories(ProMP_LearnEx
PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${MP_INCLUDE_DIRS}
)

target_link_libraries(ProMP_LearnEx
Eigen3::Eigen
)


add_executable(ProMP_LoadEx ProMP_LoadEx.cpp
ProMP.cpp
)

target_include_directories(ProMP_LoadEx
PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${MP_INCLUDE_DIRS}
)

target_link_libraries(ProMP_LoadEx
Eigen3::Eigen
)

add_library(promp ProMP.cpp)

install(TARGETS promp DESTINATION /usr/lib)