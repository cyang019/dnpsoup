#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

#EXTRA_FLAGS="-DDNPSOUP_VERBOSE=1"
EXTRA_FLAGS=""

#BUILD_TYPE=Debug
BUILD_TYPE=Release
#BUILD_TYPE=RelWithDebInfo
CC=clang++
C=clang
#CC=c++
#C=gcc


cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_C_COMPILER=$C -DCMAKE_CXX_COMPILER=$CC -GNinja .. $EXTRA_FLAGS

ninja

