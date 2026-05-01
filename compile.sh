#!/bin/bash

rm -rf build_noibm build_ibm *.vtk

cmake -S . -B build_noibm -DUSE_IBM=OFF
cmake --build build_noibm

cmake -S . -B build_ibm -DUSE_IBM=ON
cmake --build build_ibm
