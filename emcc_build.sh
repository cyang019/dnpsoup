#!/usr/bin/env bash

FILE=dist
if [ -d "$FILE" ]; then
  rm -rf $FILE
fi
TEMP_FILE=dnpsoup_cli/dist
if [ -d "$TEMP_FILE" ]; then
  rm -rf $TEMP_FILE
fi
rm dnpsoup_cli/*.a
rm -rf build
mkdir build
cd build
cmake -DUSE_EMCC=1 -DCMAKE_BUILD_TYPE=Debug -GNinja .. -DCMAKE_CXX_COMPILER=em++
ninja

# static libs
cp dnpsoup_impl/libdnpsoup_core.a ../dnpsoup_cli/
cp matrix/matrix_impl/libmatrix.a ../dnpsoup_cli/
cp ../matrix/matrix_impl/libopenblasp-r0.3.7.a ../dnpsoup_cli/

#cd ../dnpsoup_cli
#
#em++ -o dnpsoup_cli.html -O3 -s WASM=1 -std=c++1z main.cpp engine.cpp engine_impl.cpp libdnpsoup_core.a libmatrix.a libopenblasp-r0.3.7.a -I../dnpsoup_impl/include -I../matrix/matrix_impl/include -I../build
#
#mkdir dist
#mv dnpsoup_cli.html dist
#mv dnpsoup_cli.wasm dist
#mv dnpsoup_cli.js dist
#cd ..
#mv dnpsoup_cli/dist dist

# emconfigure cmake .. -D USE_EMCC=1 -DEMSCRIPTEN_GENERATE_BITCODE_STATIC_LIBRARIES=1 -DCMAKE_BUILD_TYPE=Debug -DTARGET=HASWELL -DBINARY=64
# emmake make VERBOSE=1

#emcc -O2 dnpsoup_exec.bc -o dnpsoup_exec.js

