#!/usr/bin/bash
gfortran -c src/sub/cwrapper_interface.f90 -J./obj_f -o obj_f/cwrapper_interface.o
gfortran -c -I./obj_f src/sub/input.f90 -J./obj_f -o obj_f/input.o
gfortran -I./obj_f -c src/main/GSI_mutation++.f90 -o obj_f/GSI_mutation++.o
gfortran obj_f/GSI_mutation++.o obj_f/input.o obj_f/cwrapper_interface.o -L$MPP_DIRECTORY/install/lib -lmutation++_fortran  -lmutation++ -lstdc++ -o gsi_fortran