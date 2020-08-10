if(NOT USE_EMCC)
  add_library(dnpsoup_core SHARED "")
  message(STATUS "dnpsoup_core shared library...")
else()
  add_library(dnpsoup_core STATIC "")
  message(STATUS "dnpsoup_core static library...")
endif()


if(NOT DEFINED USE_DEBUG_MODEL)
  set(USE_DEBUG_MODEL ON)
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_DEBUG)                                       
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")                  
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)                                     
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")    
endif()

#if(NOT USE_EMCC)
#  target_compile_options(dnpsoup_core PUBLIC -Wall -Wextra -Wpedantic -Werror -emit-llvm --target=wasm32-unknown-unknown-elf -S)
#else()
#  target_compile_options(dnpsoup_core PUBLIC -Wall -Wextra -Wpedantic -Werror)  
#endif()
target_compile_options(dnpsoup_core PUBLIC -Wall -Wextra -Wpedantic)

target_link_libraries(dnpsoup_core PUBLIC matrix)
message(STATUS "dnpsoup_core linked to matrix...")

if (UNIX AND NOT APPLE)
  message(STATUS "link pthread...")
  target_link_libraries(dnpsoup_core PUBLIC pthread)
endif ()

target_include_directories(dnpsoup_core
  PRIVATE
  "${PROJECT_SOURCE_DIR}/matrix/matrix_impl/include"
  PUBLIC
  "${PROJECT_SOURCE_DIR}/dnpsoup_impl/include"
  )
#target_include_directories(dnpsoup_core
#  INTERFACE
#  "${PROJECT_SOURCE_DIR}/matrix/matrix_impl/include"
#  "${PROJECT_SOURCE_DIR}/dnpsoup_impl/include"
#  )

include("${PROJECT_SOURCE_DIR}/dnpsoup_impl/include/CMakeLists.txt")
include("${PROJECT_SOURCE_DIR}/dnpsoup_impl/src/CMakeLists.txt")

