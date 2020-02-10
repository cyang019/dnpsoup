#!/usr/bin/env bash

rm -rf build
mkdir build
cd build

BUILD_TYPE=Debug
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -GNinja ..

ninja

