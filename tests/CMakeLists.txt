file(GLOB SRCS *.cpp)

# message("tests to run: ${SRCS}")

# list(FILTER SRCS EXCLUDE REGEX ".*test_dnp_runner.cpp$")
# list(FILTER SRCS EXCLUDE REGEX ".*test_pulseseq.cpp$")
# list(FILTER SRCS EXCLUDE REGEX ".*test_BDPA.cpp$")
# list(FILTER SRCS EXCLUDE REGEX ".*test_spinsys.cpp$")

add_executable(test_dnpsoup_core ${SRCS})

set(USE_COVERAGE OFF)

target_link_libraries(
  test_dnpsoup_core
  dnpsoup_core
  gtest_main
  )

# if(${CMAKE_BUILD_TYPE} EQUAL Release)
#   set(USE_COVERAGE OFF)
# endif()

if(USE_COVERAGE)
  include(cmake/CodeCoverage.cmake)

  if(APPLE)
    message("APPLE coverage...")
    target_compile_options(test_dnpsoup_core PUBLIC
      -g -O0 --coverage)
  else()
    target_compile_options(test_dnpsoup_core PUBLIC
      -Wall -fprofile-arcs -ftest-coverage
    )
    target_link_libraries(test_dnpsoup_core gcov)
  endif()
  add_custom_command(TARGET test_dnpsoup_core POST_BUILD
    COMMAND make test
    COMMAND lcov -c -d tests/CMakeFiles -o cov.info
    COMMAND genhtml cov.info -o out
    COMMAND open out/index.html
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)

else()
  target_compile_options(test_dnpsoup_core PUBLIC
    -Wall)
endif()


target_include_directories(
  test_dnpsoup_core
  PUBLIC
    ${PROJECT_BINARY_DIR}/googletest/googletest-src/googletest/include
    ${PROJECT_SOURCE_DIR}/dnpsoup_impl/include
  )

#gtest_add_tests(TARGET      test_matrix
#                TEST_SUFFIX .noArgs
#                TEST_LIST   noArgsTests
#                )
#set_tests_properties(${noArgsTests} PROPERTIES TIMEOUT 10)
add_test(NAME test_dnpsoup_core
         COMMAND test_dnpsoup_core)


