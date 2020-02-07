#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake -DUSE_EMCC=1 -DCMAKE_BUILD_TYPE=Debug -GNinja ..
ninja

cp dnpsoup_impl/libdnpsoup_core.a ../dnpsoup_cli/
cp matrix/matrix_impl/libmatrix.a ../dnpsoup_cli/
cp ../matrix/matrix_impl/libopenblasp-r0.3.7.a ../dnpsoup_cli/

#cmake -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles" -DCMAKE_TOOLCHAIN_FILE=../emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake ..
#cmake -D USE_EMCC=BOOL:ON -DCMAKE_CXX_COMPILER=emcc -DCMAKE_BUILD_TYPE=Debug -GNinja -DCMAKE_TOOLCHAIN_FILE=../../emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake ..

cd ../dnpsoup_cli

emcc -o dnpsoup_cli.html -O3 -s WASM=1 -std=c++1z engine.cpp libdnpsoup_core.a libmatrix.a libopenblasp-r0.3.7.a -I../dnpsoup_impl/include -I../matrix/matrix_impl/include -I../build
# rm -rf build
# mkdir build
# cd build
# 
# emconfigure cmake .. -D USE_EMCC=1 -DEMSCRIPTEN_GENERATE_BITCODE_STATIC_LIBRARIES=1 -DCMAKE_BUILD_TYPE=Debug -DTARGET=HASWELL -DBINARY=64
# emmake make VERBOSE=1

#emcc -O2 dnpsoup_exec.bc -o dnpsoup_exec.js

