cmake_minimum_required(VERSION 3.21)
project(IME824)

set(CMAKE_CXX_STANDARD 14)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(IME824 main.cpp LinAlg.h test.h)
target_link_libraries(IME824 ${ARMADILLO_LIBRARIES})
