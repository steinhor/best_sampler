cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(samplertest CXX)
project(samplertest_hydro CXX)
project(viscoustest CXX)
project(bulktest CXX)

add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ../software)

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
message("Libraries: msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}" )

set(SOURCE_FILES ./samplertest.cc)
SET_SOURCE_FILES_PROPERTIES(
      samplertest.cc PROPERTIES OBJECT_DEPENDS
    ../software/lib/libmsu_sampler.a
 )
 
set(SOURCE_FILES ./samplertest_hydro.cc)
SET_SOURCE_FILES_PROPERTIES(
      samplertest_hydro.cc PROPERTIES OBJECT_DEPENDS
    ../software/lib/libmsu_sampler.a
 )
 
set(SOURCE_FILES ./bulktest.cc)
 SET_SOURCE_FILES_PROPERTIES(
	bulktest.cc PROPERTIES OBJECT_DEPENDS ../software/lib/libmsu_sampler.a
)
set(SOURCE_FILES ./viscoustest.cc)
SET_SOURCE_FILES_PROPERTIES(
	viscoustest.cc PROPERTIES OBJECT_DEPENDS ../software/lib/libmsu_sampler.a
)

add_executable(samplertest ${SOURCE_FILES})

add_executable(samplertest_hydro ${SOURCE_FILES})

add_executable(viscoustest ${SOURCE_FILES})

add_executable(bulktest ${SOURCE_FILES})

add_dependencies(samplertest extern_msu_sampler)
add_dependencies(samplertest_hydro extern_msu_sampler)
add_dependencies(viscoustest extern_msu_sampler)
add_dependencies(bulktest extern_msu_sampler)
find_library(${CMAKE_CURRENT_SOURCE_DIR}/../software/lib msu_sampler)
target_link_libraries(samplertest msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(samplertest_hydro msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(viscoustest msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(bulktest msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
