#!/bin/sh

rm -rf ./build

mkdir build

cd build

cmake .. -GNinja -DENABLE_CLANG_TIDY=ON

ninja

./filter_test