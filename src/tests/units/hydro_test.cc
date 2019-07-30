#include <iostream>
#include "unittest.h"
#include "alltests.h"
#include "hydro.h"

void hydroTest(grid &gridData, parser &inputData, parallel &mpiData) {
    // ERROR IN COMPUTED SOLUTION FROM POISSON SOLVER
    double errorVal, errorTolerance;

    // POISSON SOLVER INSTANCE
    hydro *nseSolver;

    // CREATE NEW INSTANCE OF THE POISSON SOLVER USING THE GRID
#ifndef PLANAR
    nseSolver = new hydro_d3(gridData, inputData, mpiData);
#else
    nseSolver = new hydro_d2(gridData, inputData, mpiData);
#endif

    if (rootRank == 0) {
        std::cout << "\033[35mTesting NSE solver - Hydro\033[0m" << std::endl;
    }

    errorTolerance = 1.0e-10;
    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting periodic data transfer operations\033[0m";
    }

    errorVal = nseSolver->testPeriodic();
    printResult(errorVal, errorTolerance);
}
