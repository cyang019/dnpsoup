#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

#EXTRA_FLAGS="-DDNPSOUP_VERBOSE=1"
EXPTRA_FLAGS=""

BUILD_TYPE=Release
CC=clang++
C=clang
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_C_COMPILER=$C -DCMAKE_CXX_COMPILER=$CC -GNinja .. $EXTRA_FLAGS

ninja

