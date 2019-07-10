#add_library(dnpsoup_core INTERFACE)
add_library(dnpsoup_core SHARED "")

message(STATUS "dnpsoup_core shared library...")

if(NOT DEFINED USE_DEBUG_MODEL)
  set(USE_DEBUG_MODEL ON)
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_DEBUG)                                       
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")                  
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)                                     
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")    
endif()

target_compile_options(dnpsoup_core PUBLIC -Wall -Wextra -Wpedantic -Werror)  
#target_compile_options(dnpsoup_core INTERFACE -Wall -Wextra -Wpedantic -Werror)  

target_link_libraries(dnpsoup_core PUBLIC matrix)
#target_link_libraries(dnpsoup_core INTERFACE matrix)
message(STATUS "dnpsoup_core linked to matrix...")

target_include_directories(dnpsoup_core
  PRIVATE
  "${PROJECT_SOURCE_DIR}/matrix/matrix_impl/include"
  "${PROJECT_SOURCE_DIR}/dnpsoup_impl/include"
  PUBLIC
  "${PROJECT_SOURCE_DIR}/matrix/matrix_impl/include"
  "${PROJECT_SOURCE_DIR}/dnpsoup_impl/include"
  )
#target_include_directories(dnpsoup_core
#  INTERFACE
#  "${PROJECT_SOURCE_DIR}/matrix/matrix_impl/include"
#  "${PROJECT_SOURCE_DIR}/dnpsoup_impl/include"
#  )

include("${PROJECT_SOURCE_DIR}/dnpsoup_impl/include/CMakeLists.txt")
include("${PROJECT_SOURCE_DIR}/dnpsoup_impl/src/CMakeLists.txt")
