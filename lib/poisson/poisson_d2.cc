#include "poisson.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the multigrid_d2 class derived from the poisson class
 *
 ********************************************************************************************************************************************
 */
multigrid_d2::multigrid_d2(const grid &mesh, const parser &solParam): poisson(mesh, solParam) {
    blitz::TinyVector<int, 3> loBound, upBound;

    /*****************************************************************************/

    localSizeIndex = blitz::TinyVector<int, 3>(mesh.sizeIndex(0) - int(log2(inputParams.npX)),
                                               mesh.sizeIndex(1),
                                               mesh.sizeIndex(2));

    /*****************************************************************************/

    stagFull.resize(inputParams.vcDepth + 1);
    stagCore.resize(inputParams.vcDepth + 1);

    xEnd.resize(inputParams.vcDepth + 1);
    zEnd.resize(inputParams.vcDepth + 1);

    for (int i=0; i<=inputParams.vcDepth; i++) {
        loBound = 0, 0, 0;
        upBound = mgSizeArray(localSizeIndex(0) - i) - 1, 0, mgSizeArray(localSizeIndex(2) - i) - 1;
        stagCore(i) = blitz::RectDomain<3>(loBound, upBound);

        loBound = -1, -1, -1;
        upBound = stagCore(i).ubound() - loBound;
        stagFull(i) = blitz::RectDomain<3>(loBound, upBound);

        xEnd(i) = stagCore(i).ubound(0);
        zEnd(i) = stagCore(i).ubound(2);
    }

    /*****************************************************************************/

    hx.resize(inputParams.vcDepth + 1);
    hz.resize(inputParams.vcDepth + 1);

    hx2.resize(inputParams.vcDepth + 1);
    hz2.resize(inputParams.vcDepth + 1);

    hzhx.resize(inputParams.vcDepth + 1);

    for(int i=0; i<=inputParams.vcDepth; i++) {
        hx(i) = strideValues(i)*mesh.dXi;
        hz(i) = strideValues(i)*mesh.dZt;

        hx2(i) = pow(strideValues(i)*mesh.dXi, 2.0);
        hz2(i) = pow(strideValues(i)*mesh.dZt, 2.0);

        hzhx(i) = pow(strideValues(i), 4.0)*pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);
    }

    /*****************************************************************************/

    pressureData.resize(inputParams.vcDepth + 1);
    residualData.resize(inputParams.vcDepth + 1);
    tmpDataArray.resize(inputParams.vcDepth + 1);
    smoothedPres.resize(inputParams.vcDepth + 1);

    for (int i=0; i <= inputParams.vcDepth; i++) {
        pressureData(i).resize(blitz::TinyVector<int, 3>(stagFull(i).ubound() - stagFull(i).lbound() + 1));
        pressureData(i).reindexSelf(stagFull(i).lbound());
        pressureData(i) = 0.0;

        tmpDataArray(i).resize(blitz::TinyVector<int, 3>(stagFull(i).ubound() - stagFull(i).lbound() + 1));
        tmpDataArray(i).reindexSelf(stagFull(i).lbound());
        tmpDataArray(i) = 0.0;

        residualData(i).resize(blitz::TinyVector<int, 3>(stagFull(i).ubound() - stagFull(i).lbound() + 1));
        residualData(i).reindexSelf(stagFull(i).lbound());
        residualData(i) = 0.0;

        smoothedPres(i).resize(blitz::TinyVector<int, 3>(stagFull(i).ubound() - stagFull(i).lbound() + 1));
        smoothedPres(i).reindexSelf(stagFull(i).lbound());
        smoothedPres(i) = 0.0;
    }

    initDirichlet();
}


void multigrid_d2::computeResidual() {
    tmpDataArray(vLevel) = 0.0;

    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            tmpDataArray(vLevel)(i, 0, k) = residualData(vLevel)(i, 0, k) -
                                          ((pressureData(vLevel)(i + 1, 0, k) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i - 1, 0, k))/(hx2(vLevel)) +
                                           (pressureData(vLevel)(i, 0, k + 1) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i, 0, k - 1))/(hz2(vLevel)));
        }
    }
}


void multigrid_d2::solve() {
    int iterCount = 0;

    while (true) {
        imposeBC();

        // GAUSS-SEIDEL ITERATIVE SOLVER
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                pressureData(vLevel)(i, 0, k) = (hz2(vLevel) * (pressureData(vLevel)(i + 1, 0, k) + pressureData(vLevel)(i - 1, 0, k)) +
                                                 hx2(vLevel) * (pressureData(vLevel)(i, 0, k + 1) + pressureData(vLevel)(i, 0, k - 1)) -
                                                hzhx(vLevel) *  residualData(vLevel)(i, 0, k))/
                                            (2.0 * (hz2(vLevel) + hx2(vLevel)));
            }
        }

        real tempValue = 0.0;
        real globalMax = -1.0e-10;
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                tempValue = fabs(residualData(vLevel)(i, 0, k) -
                               ((pressureData(vLevel)(i + 1, 0, k) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i - 1, 0, k))/(hx2(vLevel)) +
                                (pressureData(vLevel)(i, 0, k + 1) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i, 0, k - 1))/(hz2(vLevel))));

                if (tempValue > globalMax) {
                    globalMax = tempValue;
                }
            }
        }

        std::cout << vLevel << "\t" << iterCount << "\t" << globalMax << std::endl;

        if (globalMax < 1.0e-6) {
            break;
        }

        iterCount += 1;
        if (iterCount > maxCount) {
            std::cout << "ERROR: Jacobi iterations for solution at coarsest level not converging. Aborting" << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }

    imposeBC();
}


void multigrid_d2::smooth(const int smoothCount) {
    for(int n=0; n<smoothCount; n++) {
        imposeBC();

        // GAUSS-SEIDEL ITERATIVE SMOOTHING
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                pressureData(vLevel)(i, 0, k) = (hz2(vLevel) * (pressureData(vLevel)(i + 1, 0, k) + pressureData(vLevel)(i - 1, 0, k)) +
                                                 hx2(vLevel) * (pressureData(vLevel)(i, 0, k + 1) + pressureData(vLevel)(i, 0, k - 1)) -
                                                hzhx(vLevel) *  residualData(vLevel)(i, 0, k))/
                                            (2.0 * (hz2(vLevel) + hx2(vLevel)));
            }
        }
    }

    imposeBC();
}


void multigrid_d2::coarsen() {
    real facePoints, vertPoints;

    int i2, k2;
    int pLevel;

    pLevel = vLevel;
    vLevel += 1;

    residualData(vLevel) = 0.0;

    for (int i = 0; i <= xEnd(vLevel); ++i) {
        i2 = i*2;
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            k2 = k*2;
            facePoints = (tmpDataArray(pLevel)(i2 + 1, 0, k2) + tmpDataArray(pLevel)(i2 - 1, 0, k2) +
                          tmpDataArray(pLevel)(i2, 0, k2 + 1) + tmpDataArray(pLevel)(i2, 0, k2 - 1))*0.125;
            vertPoints = (tmpDataArray(pLevel)(i2 + 1, 0, k2 + 1) +
                          tmpDataArray(pLevel)(i2 + 1, 0, k2 - 1) +
                          tmpDataArray(pLevel)(i2 - 1, 0, k2 + 1) +
                          tmpDataArray(pLevel)(i2 - 1, 0, k2 - 1))*0.0625;

            residualData(vLevel)(i, 0, k) = facePoints + vertPoints + tmpDataArray(pLevel)(i2, 0, k2)*0.25;
        }
    }
}


void multigrid_d2::prolong() {
    int i2, k2;
    int pLevel;

    pLevel = vLevel;
    vLevel -= 1;

    pressureData(vLevel) = 0.0;

    for (int i = 0; i <= xEnd(vLevel); i++) {
        i2 = i/2;
        if (isOdd(i)) {
            for (int k = 0; k <= zEnd(vLevel); k++) {
                k2 = k/2;
                if (isOdd(k)) { // Both i and k are odd
                    pressureData(vLevel)(i, 0, k) = (pressureData(pLevel)(i2, 0, k2)     + pressureData(pLevel)(i2, 0, k2 + 1) +
                                                     pressureData(pLevel)(i2 + 1, 0, k2) + pressureData(pLevel)(i2 + 1, 0, k2 + 1))/4.0;
                } else {        // Here i is odd, but k is even
                    pressureData(vLevel)(i, 0, k) = (pressureData(pLevel)(i2, 0, k2) + pressureData(pLevel)(i2 + 1, 0, k2))/2.0;
                }
            }
        } else {
            for (int k = 0; k <= zEnd(vLevel); k++) {
                k2 = k/2;
                if (isOdd(k)) { // Here i is even, but k is odd
                    pressureData(vLevel)(i, 0, k) = (pressureData(pLevel)(i2, 0, k2) + pressureData(pLevel)(i2, 0, k2 + 1))/2.0;
                } else {        // Both i and k are even
                    pressureData(vLevel)(i, 0, k) = pressureData(pLevel)(i2, 0, k2);
                }
            }
        }
    }
}


real multigrid_d2::computeError(const int normOrder) {
    real residualVal = 0.0;

    real tempValue = 0.0;
    real numValLoc = 0.0;
    real denValLoc = 0.0;
    int valCountLoc = 0;

    for (int i = 0; i <= xEnd(0); ++i) {
        for (int k = 0; k <= zEnd(0); ++k) {
            tempValue = fabs(((pressureData(0)(i + 1, 0, k) - 2.0*pressureData(0)(i, 0, k) + pressureData(0)(i - 1, 0, k))/hx2(0) +
                              (pressureData(0)(i, 0, k + 1) - 2.0*pressureData(0)(i, 0, k) + pressureData(0)(i, 0, k - 1))/hz2(0)) -
                               residualData(0)(i, 0, k));

            if (normOrder == 0) {
                if (tempValue > numValLoc) numValLoc = tempValue;
            } else {
                numValLoc += tempValue*tempValue;
                denValLoc += residualData(0)(i, 0, k)*residualData(0)(i, 0, k);
                valCountLoc += 1;
            }
        }
    }

    if (normOrder == 0) {
        denValLoc = blitz::max(fabs(residualData(0)));
        residualVal = numValLoc/denValLoc;
        //residualVal = numValLoc;
    } else {
        residualVal = sqrt(numValLoc/valCountLoc)/sqrt(denValLoc/valCountLoc);
        //residualVal = sqrt(numValLoc/valCountLoc);
    }

    return residualVal;
}


void multigrid_d2::initDirichlet() {
    real xDist, zDist;

    xWall.resize(stagFull(0).ubound(2) - stagFull(0).lbound(2) + 1);
    xWall.reindexSelf(stagFull(0).lbound(2));
    xWall = 0.0;

    zWall.resize(stagFull(0).ubound(0) - stagFull(0).lbound(0) + 1);
    zWall.reindexSelf(stagFull(0).lbound(0));
    zWall = 0.0;

    // Compute values at the walls using the (r^2)/4 formula
    xDist = hx(0)*(int(mgSizeArray(localSizeIndex(0))/2));
    for (int k=stagCore(0).lbound(2); k<=stagCore(0).ubound(2); k++) {
        zDist = hz(0)*(k - stagCore(0).ubound(2)/2);
        xWall(k) = (xDist*xDist + zDist*zDist)/4.0;
    }

    zDist = hz(0)*(int(mgSizeArray(localSizeIndex(2))/2));
    for (int i=stagCore(0).lbound(0); i<=stagCore(0).ubound(0); i++) {
        xDist = hx(0)*(i - stagCore(0).ubound(0)/2);
        zWall(i) = (xDist*xDist + zDist*zDist)/4.0;
    }
}


void multigrid_d2::imposeBC() {
    if (zeroBC) {
        pressureData(vLevel)(-1, 0, all) = -pressureData(vLevel)(1, 0, all);

        pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, 0, all) = -pressureData(vLevel)(stagCore(vLevel).ubound(0) - 1, 0, all);
    } else {
        pressureData(vLevel)(-1, 0, all) = xWall(all);

        pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, 0, all) = xWall(all);
    }

    if (zeroBC) {
        pressureData(vLevel)(all, 0, -1) = -pressureData(vLevel)(all, 0, 1);

        pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) + 1) = -pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) - 1);
    } else {
        pressureData(vLevel)(all, 0, -1) = zWall(all);

        pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) + 1) = zWall(all);
    }
}
