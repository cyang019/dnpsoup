cmake_minimum_required (VERSION 3.14)
project (dnpsoup 
  LANGUAGES CXX
  VERSION 0.1)

message(STATUS "Project dnpsoup...")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# if(NOT DEFINED CMAKE_CXX_STANDARD)
#   set(CMAKE_CXX_STANDARD 17)
#   message(STATUS "[dnpsoup] cpp standard - ${CMAKE_CXX_STANDARD}")
# endif()

# set(CMAKE_CXX_COMPILER /usr/bin/clang++)
# message(STATUS "[dnpsoup] compiler - ${CMAKE_CXX_COMPILER}")

if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/configure_dnpsoup.h.in ${PROJECT_BINARY_DIR}/configure_dnpsoup.h)
message( STATUS "configure_dnpsoup.h => ${PROJECT_BINARY_DIR}" )
include_directories(${PROJECT_BINARY_DIR})

add_subdirectory("matrix/matrix_impl")
add_subdirectory("dnpsoup_impl")
add_subdirectory("dnpsoup_cli")
option(BUILD_TESTS "build unittests as well for dnpsoup" ON)
if(BUILD_TESTS)
  enable_testing()
  message(STATUS "build unittests for dnpsoup...")
  add_subdirectory(tests)

  include(cmake/googletest.cmake)
  fetch_googletest(
    ${PROJECT_SOURCE_DIR}/cmake
    ${PROJECT_BINARY_DIR}/googletest
    )
endif()

