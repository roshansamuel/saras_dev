#include <iostream>
#include "unittest.h"
#include "alltests.h"
#include "poisson.h"
#include "sfield.h"

static void taylorGreen(sfield &rho, sfield &P_analytic, grid &mesh);

void poissonTest(grid &gridData, parser &inputData) {
    // ERROR IN COMPUTED SOLUTION FROM POISSON SOLVER
    double errorVal, errorTolerance;

    // POISSON SOLVER INSTANCE
    poisson *mgSolver;

    // SCALAR FIELDS TO TEST POISSON SOLVER
    sfield rho(gridData, "rho", true, true, true);
    sfield P_analytic(gridData, "P_anl", true, true, true);
    sfield P_calculat(gridData, "P_cal", true, true, true);

    // CREATE NEW INSTANCE OF THE POISSON SOLVER USING THE GRID
#ifndef PLANAR
    mgSolver = new multigrid_d3(gridData, inputData);
#else
    mgSolver = new multigrid_d2(gridData, inputData);
#endif

    if (rootRank == 0) {
        std::cout << "\033[35mTesting Poisson solver\033[0m" << std::endl;
    }

    errorTolerance = 1.0e-10;
    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting prolongation operations\033[0m";
    }

    errorVal = mgSolver->testProlong();
    printResult(errorVal, errorTolerance);

    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting strided data transfer operations\033[0m";
    }

    errorVal = mgSolver->testTransfer();
    printResult(errorVal, errorTolerance);

    taylorGreen(rho, P_analytic, gridData);

    if (inputData.probType == 2) {
        if (rootRank == 0) {
            std::cout << std::setw(60) << std::left << "\033[34mTesting periodic data transfer operations\033[0m";
        }

        errorVal = mgSolver->testPeriodic();
        printResult(errorVal, errorTolerance);
    }
}

static void taylorGreen(sfield &rho, sfield &P_analytic, grid &mesh) {
    // ANALYTIC SOLUTION AND RHS FOR THE POISSON SOLVER
#ifndef PLANAR
    for (int i=rho.F.F.lbound(0); i <= rho.F.F.ubound(0); i++) {
        for (int j=rho.F.F.lbound(1); j <= rho.F.F.ubound(1); j++) {
            for (int k=rho.F.F.lbound(2); k <= rho.F.F.ubound(2); k++) {
                P_analytic.F.F(i, j, k) = sin(1.0*M_PI*mesh.xStaggr(i))*
                                          cos(2.0*M_PI*mesh.yStaggr(j))*
                                          cos(4.0*M_PI*mesh.zStaggr(k));

                rho.F.F(i, j, k) = -21.0*M_PI*M_PI*sin(1.0*M_PI*mesh.xStaggr(i))*
                                                   cos(2.0*M_PI*mesh.yStaggr(j))*
                                                   cos(4.0*M_PI*mesh.zStaggr(k));
            }
        }
    }
#else
    int j = 0;
    for (int i=rho.F.F.lbound(0); i <= rho.F.F.ubound(0); i++) {
        for (int k=rho.F.F.lbound(2); k <= rho.F.F.ubound(2); k++) {
            P_analytic.F.F(i, j, k) = sin(1.0*M_PI*mesh.xStaggr(i))*
                                      cos(4.0*M_PI*mesh.zStaggr(k));

            rho.F.F(i, j, k) = -17.0*M_PI*M_PI*sin(1.0*M_PI*mesh.xStaggr(i))*
                                               cos(4.0*M_PI*mesh.zStaggr(k));
        }
    }
#endif
}
