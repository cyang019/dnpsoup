#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

#cmake -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles" -DCMAKE_TOOLCHAIN_FILE=../emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake ..
#cmake -D USE_EMCC=BOOL:ON -DCMAKE_CXX_COMPILER=emcc -DCMAKE_BUILD_TYPE=Debug -GNinja -DCMAKE_TOOLCHAIN_FILE=../../emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake ..
cmake -D USE_EMCC=BOOL:ON -DCMAKE_BUILD_TYPE=Debug -GNinja ..
ninja

cp dnpsoup_impl/libdnpsoup_core.a ../dnpsoup_cli/
cd ../dnpsoup_cli
rm -rf build
mkdir build
cd build

emconfigure cmake .. -D USE_EMCC=BOOL:ON
emmake make VERBOSE=1

emcc -O2 dnpsoup_exec.bc -o dnpsoup_exec.js

