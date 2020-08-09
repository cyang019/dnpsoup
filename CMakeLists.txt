cmake_minimum_required (VERSION 3.0.2)
project (dnpsoup 
  LANGUAGES CXX
  VERSION 0.1)
enable_language(CXX)

message(STATUS "Project dnpsoup...")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

configure_file(${CMAKE_CURRENT_LIST_DIR}/configure_dnpsoup.h.in ${PROJECT_BINARY_DIR}/configure_dnpsoup.h)
message( STATUS "configure_dnpsoup.h => ${PROJECT_BINARY_DIR}" )
include_directories(${PROJECT_BINARY_DIR})

#if(USE_EMCC)
  #message(STATUS "build OpenBLAS from source...")
  #set(OpenBLAS_VERSION "0.3.7")
  #add_subdirectory("OpenBLAS-${OpenBLAS_VERSION}")
#endif()
add_subdirectory("matrix/matrix_impl")
add_subdirectory("dnpsoup_impl")

set(CMAKE_VERBOSE_MAKEFILE ON)
if(NOT USE_EMCC)
  add_subdirectory("dnpsoup_cli")
  option(BUILD_TESTS "build unittests as well for dnpsoup" ON)
else()
  option(BUILD_TESTS "build unittests as well for dnpsoup" OFF)
endif()

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

