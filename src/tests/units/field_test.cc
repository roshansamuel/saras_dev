#include <iostream>
#include "unittest.h"
#include "alltests.h"
#include "sfield.h"
#include "vfield.h"

static void taylorGreen(vfield &V, grid &mesh);

void fieldTest(grid &gridData) {
    // ERROR IN COMPUTED DIVERGENCE
    double errorVal, errorTolerance;

    sfield div(gridData, "DIV", true, true, true);
    vfield V(gridData, "V");

    // TOLERANCE IN ERROR OF COMPUTED VALUES
#ifdef PLANAR
    errorTolerance = std::max(gridData.dXi, gridData.dZt);
#else
    errorTolerance = std::max(std::max(gridData.dXi, gridData.dEt), gridData.dZt);
#endif

    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting divergence calculation of vfield\033[0m";
    }

    taylorGreen(V, gridData);

    V.divergence(div);

    errorVal = div.fieldMax();

    printResult(errorVal, errorTolerance);
}

static void taylorGreen(vfield &V, grid &mesh) {
    // X-Velocity
#ifndef PLANAR
    for (int i=V.Vx.F.F.lbound(0); i <= V.Vx.F.F.ubound(0); i++) {
        for (int j=V.Vx.F.F.lbound(1); j <= V.Vx.F.F.ubound(1); j++) {
            for (int k=V.Vx.F.F.lbound(2); k <= V.Vx.F.F.ubound(2); k++) {
                V.Vx.F.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                    cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                    cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }
#else
    int j = 0;
    for (int i=V.Vx.F.F.lbound(0); i <= V.Vx.F.F.ubound(0); i++) {
        for (int k=V.Vx.F.F.lbound(2); k <= V.Vx.F.F.ubound(2); k++) {
            V.Vx.F.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
        }
    }
#endif

    // Y-Velocity
#ifndef PLANAR
    for (int i=V.Vy.F.F.lbound(0); i <= V.Vy.F.F.ubound(0); i++) {
        for (int j=V.Vy.F.F.lbound(1); j <= V.Vy.F.F.ubound(1); j++) {
            for (int k=V.Vy.F.F.lbound(2); k <= V.Vy.F.F.ubound(2); k++) {
                V.Vy.F.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                     sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                     cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }
#else
    V.Vy.F.F = 0.0;
#endif

    // Z-Velocity
#ifndef PLANAR
    V.Vz.F.F = 0.0;
#else
    for (int i=V.Vz.F.F.lbound(0); i <= V.Vz.F.F.ubound(0); i++) {
        for (int k=V.Vz.F.F.lbound(2); k <= V.Vz.F.F.ubound(2); k++) {
            V.Vz.F.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                 sin(2.0*M_PI*mesh.zColloc(k)/mesh.zLen);
        }
    }
#endif
}
