cmake_minimum_required (VERSION 3.6.0)
project (dnpsoup 
  LANGUAGES CXX
  VERSION 0.1.0)
enable_language(CXX)

message(STATUS "Project dnpsoup...")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(MSVC)
  message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
  if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/out")
  endif()
  if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/out")
  endif()
  if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/out")
  endif()
endif()

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

  include(cmake/googletest.cmake)
  fetch_googletest(
    ${PROJECT_SOURCE_DIR}/cmake
    ${PROJECT_BINARY_DIR}/googletest
    )

  include_directories(${GTEST_INCLUDE_DIRS})

  add_subdirectory(tests)
endif()

#if (BUILD_TESTS)
#  install(TARGETS matrix dnpsoup_core dnpsoup_exec test_dnpsoup_core
#          RUNTIME DESTINATION bin
#          LIBRARY DESTINATION lib
#          ARCHIVE DESTINATION lib)
#else()
#  install(TARGETS matrix dnpsoup_core dnpsoup_exec
#          RUNTIME DESTINATION bin
#          LIBRARY DESTINATION lib
#          ARCHIVE DESTINATION lib)
#endif()

