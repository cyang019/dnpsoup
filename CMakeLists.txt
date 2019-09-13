cmake_minimum_required (VERSION 3.14)
project (dnpsoup CXX)

message(STATUS "Project dnpsoup...")

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  message(STATUS "cpp standard - ${CMAKE_CXX_STANDARD}")
endif()

if(APPLE)
  set(CMAKE_CXX_COMPILER /usr/bin/clang++)
  message(STATUS "compiler - ${CMAKE_CXX_COMPILER}")
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()


add_subdirectory("matrix/matrix_impl")
add_subdirectory("dnpsoup_impl")
option(BUILD_TESTS "build unittests as well for dnpsoup" ON)
if(BUILD_TESTS)
  enable_testing()
  message(STATUS "build unittests for dnpsoup...")
  add_subdirectory(tests)
endif()

include(cmake/googletest.cmake)
include(cmake/CodeCoverage.cmake)
fetch_googletest(
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_BINARY_DIR}/googletest
  )

