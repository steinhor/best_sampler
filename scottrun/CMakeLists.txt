cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(samplertest CXX)

add_custom_target(extern_pratt_sampler COMMAND make WORKING_DIRECTORY ../software)

#include(${CMAKE_CURRENT_SOURCE_DIR}/../software/CMakeLists.txt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../software/include")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../software/lib)
find_library(${CMAKE_CURRENT_SOURCE_DIR}/../software/lib pratt_sampler)
message("Libraries: pratt_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/samplertest.cc)

add_executable(samplertest ${SOURCE_FILES})
target_link_libraries(samplertest pratt_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
add_dependencies(samplertest extern_pratt_sampler)