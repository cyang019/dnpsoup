#add_library(dnpsoup_core INTERFACE)
add_executable(dnpsoup_exec engine.h engine.cpp main.cpp)

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

target_link_libraries(dnpsoup_exec PRIVATE dnpsoup_core)
#target_link_libraries(dnpsoup_core INTERFACE matrix)
message(STATUS "ex_dnpsoup linked to dnpsoup_core...")

target_include_directories(dnpsoup_exec
  PRIVATE
  "${PROJECT_SOURCE_DIR}/matrix/matrix_impl/include"
  "${PROJECT_SOURCE_DIR}/dnpsoup_impl/include"
  "${PROJECT_SOURCE_DIR}/dnpsoup_cli/include"
  )

set(USE_EMCC OFF)
if(USE_EMCC)
  message(STATUS "use emcc...")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist)
  set(CMAKE_BUILD_TYPE_INIT "Release")
  set(CMAKE_C_COMPILER "emcc")
  set(CMAKE_CXX_COMPILER "emcc")

  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

  include_directories(
    "${PROECT_SOURCE_DIR}/../emsdk/upstream/include"
    "${PROECT_SOURCE_DIR}/../emsdk/upstream/emscripten/system/include"
    )

  set_target_properties(dnpsoup_exec PROPERTIES LINK_FLAGS "-o dist/dnpsoup.js")
endif()

