#!/bin/bash

# Test for 2D LDC and comparison with Ghia et al's (1982, J. Comput. Phys., 48, 387 - 411) result
PROC=4

# If build directory doesn't exist, create it
if [ ! -d build ]; then
    mkdir build
fi

# Switch to build directory
cd build

# Run cmake with necessary flags for 2D LDC test
CC=mpicc CXX=mpicxx cmake ../../ -DPLANAR=ON -DREAL_DOUBLE=ON

# Compile
make -j8

# Move the executable to the directory where the test will be performed
mv ../../saras ../../tests/ldcTest/

# Switch to ldcTest directory
cd ../../tests/ldcTest/

# Run the test case
mpirun -np $PROC ./saras
#CC=mpicc CXX=mpicxx cmake ../../ -DTEST_RUN=ON -DPLANAR=ON -DREAL_SINGLE=ON

# Run the python script to read the output file and compare with Ghia results
python validate_ldc.py
