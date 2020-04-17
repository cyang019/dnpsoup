if(USE_EMCC)
  cmake_minimum_required(VERSION 3.14)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

add_executable(dnpsoup_exec engine.h engine_impl.h engine.cpp engine_impl.cpp main.cpp)

message(STATUS "dnpsoup_exec...")

if(NOT DEFINED USE_DEBUG_MODEL)
  set(USE_DEBUG_MODEL ON)
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_DEBUG)                                       
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")                  
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)                                     
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")    
endif()

target_compile_options(dnpsoup_exec PUBLIC -Wall -Wextra -Wpedantic -Werror)  
#target_compile_options(dnpsoup_core INTERFACE -Wall -Wextra -Wpedantic -Werror)  

if(USE_EMCC)
  target_link_libraries(dnpsoup_exec PRIVATE "${CMAKE_CURRENT_LIST_DIR}/dnpsoup_core")
else()
  target_link_libraries(dnpsoup_exec PRIVATE dnpsoup_core)
endif()
#target_link_libraries(dnpsoup_core INTERFACE matrix)
message(STATUS "ex_dnpsoup linked to dnpsoup_core...")

target_include_directories(dnpsoup_exec
  PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/../matrix/matrix_impl/include"
  "${CMAKE_CURRENT_LIST_DIR}/../dnpsoup_impl/include"
  "${CMAKE_CURRENT_LIST_DIR}/../dnpsoup_cli/include"
  )

if(USE_EMCC)
  include_directories("${CMAKE_CURRENT_LIST_DIR}/../build")
endif()

#if(USE_EMCC)
#  message(STATUS "use emcc...")
#  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist)
#  set(CMAKE_BUILD_TYPE_INIT "Release")
#  set(CMAKE_C_COMPILER "emcc")
#  set(CMAKE_CXX_COMPILER "emcc")
#
#  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
#
#  include_directories(
#    "${PROECT_SOURCE_DIR}/../emsdk/upstream/include"
#    "${PROECT_SOURCE_DIR}/../emsdk/upstream/emscripten/system/include"
#    )
#
#  set_target_properties(dnpsoup_exec PROPERTIES LINK_FLAGS "-o dist/dnpsoup.js")
#endif()

