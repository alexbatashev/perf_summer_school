cmake_minimum_required(VERSION 3.17)

project(
  summer_school
  VERSION 1.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

if (UNIX)
  set(CMAKE_CXX_FLAGS "-fno-omit-frame-pointer ${CMAKE_CXX_FLAGS}")
  if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "-gline-tables-only ${CMAKE_CXX_FLAGS}")
  endif()

endif()

include(cmake/CPM.cmake)

set(BENCHMARK_ENABLE_WERROR OFF)
set(BENCHMARK_ENABLE_TESTING OFF)

CPMAddPackage("gh:google/benchmark@1.8.4")
CPMAddPackage("gh:catchorg/Catch2@3.6.0")

add_library(task1_lib STATIC tasks/task1.cpp)
add_executable(task1 tasks/task1_bench.cpp)
target_link_libraries(task1 PRIVATE task1_lib benchmark::benchmark)
target_include_directories(task1_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tasks)

add_library(task2_lib STATIC tasks/task2.cpp)
add_executable(task2 tasks/task2_bench.cpp)
target_link_libraries(task2 PRIVATE task2_lib benchmark::benchmark)
target_include_directories(task2_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tasks)

add_executable(tests tests/task1.cpp tests/task2.cpp)
target_link_libraries(tests PRIVATE task1_lib task2_lib Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(tests)
