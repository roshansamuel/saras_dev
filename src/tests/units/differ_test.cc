#include <iostream>
#include "unittest.h"
#include "alltests.h"
#include "sfield.h"

static void taylorGreen(sfield &F, sfield &dF, grid &mesh);

void differTest(grid &gridData) {
    // ERROR IN COMPUTED DIVERGENCE
    double errorTolerance;

    sfield F(gridData, "F", true, true, true);
    sfield dF(gridData, "dF", true, true, true);

    errorTolerance = gridData.dXi*gridData.dXi*4.0*M_PI*M_PI;

    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting computation of derivatives in sfield\033[0m";
    }

    taylorGreen(F, dF, gridData);

    F.calcDerivatives1();

    testError(dF.F.F, F.dF_dx, 1, errorTolerance);
}

static void taylorGreen(sfield &F, sfield &dF, grid &mesh) {
#ifndef PLANAR
    for (int i=F.F.F.lbound(0); i <= F.F.F.ubound(0); i++) {
        for (int j=F.F.F.lbound(1); j <= F.F.F.ubound(1); j++) {
            for (int k=F.F.F.lbound(2); k <= F.F.F.ubound(2); k++) {
                F.F.F(i, j, k) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                 cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                 cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                dF.F.F(i, j, k) = 2.0*M_PI*cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                           cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                           cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }
#else
    int j = 0;
    for (int i=F.F.F.lbound(0); i <= F.F.F.ubound(0); i++) {
        for (int k=F.F.F.lbound(2); k <= F.F.F.ubound(2); k++) {
            F.F.F(i, j, k) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                             cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            dF.F.F(i, j, k) = 2.0*M_PI*cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                       cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
        }
    }
#endif
}
