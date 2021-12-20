#!/bin/bash
rm -rf build
mkdir build
cd build
cmake -DENABLE_TESTING=ON -DBUILD_FORTRAN_WRAPPER=ON -DBUILD_PYTHON_WRAPPER=ON -DCMAKE_INSTALL_PREFIX:PATH=$(realpath ../install) ..
make -j
make install
