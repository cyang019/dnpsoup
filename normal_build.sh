#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

#EXTRA_FLAGS="-DDNPSOUP_VERBOSE=1"
EXPTRA_FLAGS=""

BUILD_TYPE=Release
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -GNinja .. $EXTRA_FLAGS

ninja

