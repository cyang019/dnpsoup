#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Debug -GNinja -DCMAKE_TOOLCHAIN_FILE=../emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake ..

