#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

EXTRA_FLAG="-DDNPSOUP_VERBOSE=1"

BUILD_TYPE=Debug
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -GNinja .. $EXTRA_FLAG

ninja

