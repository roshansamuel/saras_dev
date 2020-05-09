/********************************************************************************************************************************************
 * Saras
 * 
 * Copyright (C) 2019, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file poisson3.cc
 *
 *  \brief Definitions for functions of class poisson for 3D
 *  \sa poisson.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "poisson.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the multigrid_d3 class derived from the poisson class
 *
 *          The constructor of the derived multigrid_d3 class frst calls the base poisson class with the arguments passed to it.
 *          It then calls a series of functions in sequence to initialize all the necessary parameters and data structures to
 *          store and manipulate the multi-grid data.
 *          Since the multi-grid solver operates on the staggered grid, it first computes the limits of the full and core
 *          staggered grid, as the grid class does the same for the collocated grid.
 *
 *          It then initializes all the Range objects to obtain the correct slices of the full grid at various
 *          levels of the V-cycle.
 *          It also copies the staggered grid derivatives to local arrays with wide pads, and finally generates the MPI datatypes
 *          for data transfer between sub-domain boundaries.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
multigrid_d3::multigrid_d3(const grid &mesh, const parser &solParam): poisson(mesh, solParam) {
    // GET THE localSizeIndex AS IT WILL BE USED TO SET THE FULL AND CORE LIMITS OF THE STAGGERED POINTS
    setLocalSizeIndex();

    // SET THE FULL AND CORE LIMTS SET ABOVE USING THE localSizeIndex VARIABLE SET ABOVE
    setStagBounds();

    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // COPY THE STAGGERED GRID DERIVATIVES TO LOCAL ARRAYS
    copyStaggrDerivs();

    // RESIZE AND INITIALIZE NECESSARY DATA-STRUCTURES
    initializeArrays();

    // INITIALIZE DIRICHLET BCs WHEN TESTING THE POISSON SOLVER
#ifdef TEST_POISSON
    initDirichlet();
#endif
}


void multigrid_d3::computeResidual() {
    // Calculate Laplacian of the pressure field and subtract it from the RHS of Poisson equation to obtain the residual
    // Needed update: Substitute the below OpenMP parallel loop with vectorized Blitz operation and check for speed increase.
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
    for (int iX = 0; iX <= xEnd(vLevel); iX += 1) {
        for (int iY = 0; iY <= yEnd(vLevel); iY += 1) {
            for (int iZ = 0; iZ <= zEnd(vLevel); iZ += 1) {
                tmpDataArray(vLevel)(iX, iY, iZ) =  residualData(vLevel)(iX, iY, iZ) -
                               (xix2(vLevel)(iX) * (pressureData(vLevel)(iX + 1, iY, iZ) - 2.0*pressureData(vLevel)(iX, iY, iZ) + pressureData(vLevel)(iX - 1, iY, iZ))/(hx(vLevel)*hx(vLevel)) +
                                xixx(vLevel)(iX) * (pressureData(vLevel)(iX + 1, iY, iZ) - pressureData(vLevel)(iX - 1, iY, iZ))/(2.0*hx(vLevel)) +
                                ety2(vLevel)(iY) * (pressureData(vLevel)(iX, iY + 1, iZ) - 2.0*pressureData(vLevel)(iX, iY, iZ) + pressureData(vLevel)(iX, iY - 1, iZ))/(hy(vLevel)*hy(vLevel)) +
                                etyy(vLevel)(iY) * (pressureData(vLevel)(iX, iY + 1, iZ) - pressureData(vLevel)(iX, iY - 1, iZ))/(2.0*hy(vLevel)) +
                                ztz2(vLevel)(iZ) * (pressureData(vLevel)(iX, iY, iZ + 1) - 2.0*pressureData(vLevel)(iX, iY, iZ) + pressureData(vLevel)(iX, iY, iZ - 1))/(hz(vLevel)*hz(vLevel)) +
                                ztzz(vLevel)(iZ) * (pressureData(vLevel)(iX, iY, iZ + 1) - pressureData(vLevel)(iX, iY, iZ - 1))/(2.0*hz(vLevel)));
            }
        }
    }
}


void multigrid_d3::smooth(const int smoothCount) {
#ifdef TIME_RUN
    struct timeval begin, end;
#endif
    tmpDataArray(vLevel) = 0.0;

    for(int n=0; n<smoothCount; n++) {
#ifdef TIME_RUN
        gettimeofday(&begin, NULL);
#endif
        // IMPOSE BOUNDARY CONDITION
        imposeBC();

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        smothTimeTran += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

        gettimeofday(&begin, NULL);
#endif

        if (inputParams.gsSmooth) {
            // GAUSS-SEIDEL ITERATIVE SMOOTHING
            for (int iX = 0; iX <= xEnd(vLevel); iX += 1) {
                for (int iY = 0; iY <= yEnd(vLevel); iY += 1) {
                    for (int iZ = 0; iZ <= zEnd(vLevel); iZ += 1) {
                        tmpDataArray(vLevel)(iX, iY, iZ) = (hyhz(vLevel) * xix2(vLevel)(iX) * (pressureData(vLevel)(iX + 1, iY, iZ) + tmpDataArray(vLevel)(iX - 1, iY, iZ))*2.0 +
                                                            hyhz(vLevel) * xixx(vLevel)(iX) * (pressureData(vLevel)(iX + 1, iY, iZ) - tmpDataArray(vLevel)(iX - 1, iY, iZ))*hx(vLevel) +
                                                            hzhx(vLevel) * ety2(vLevel)(iY) * (pressureData(vLevel)(iX, iY + 1, iZ) + tmpDataArray(vLevel)(iX, iY - 1, iZ))*2.0 +
                                                            hzhx(vLevel) * etyy(vLevel)(iY) * (pressureData(vLevel)(iX, iY + 1, iZ) - tmpDataArray(vLevel)(iX, iY - 1, iZ))*hy(vLevel) +
                                                            hxhy(vLevel) * ztz2(vLevel)(iZ) * (pressureData(vLevel)(iX, iY, iZ + 1) + tmpDataArray(vLevel)(iX, iY, iZ - 1))*2.0 +
                                                            hxhy(vLevel) * ztzz(vLevel)(iZ) * (pressureData(vLevel)(iX, iY, iZ + 1) - tmpDataArray(vLevel)(iX, iY, iZ - 1))*hz(vLevel) -
                                                    2.0 * hxhyhz(vLevel) * residualData(vLevel)(iX, iY, iZ))/
                                                   (4.0 * (hyhz(vLevel)*xix2(vLevel)(iX) + hzhx(vLevel)*ety2(vLevel)(iY) + hxhy(vLevel)*ztz2(vLevel)(iZ)));
                    }
                }
            }
        } else {
            // JACOBI ITERATIVE SMOOTHING
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
            for (int iX = 0; iX <= xEnd(vLevel); iX += 1) {
                for (int iY = 0; iY <= yEnd(vLevel); iY += 1) {
                    for (int iZ = 0; iZ <= zEnd(vLevel); iZ += 1) {
                        tmpDataArray(vLevel)(iX, iY, iZ) = (hyhz(vLevel) * xix2(vLevel)(iX) * (pressureData(vLevel)(iX + 1, iY, iZ) + pressureData(vLevel)(iX - 1, iY, iZ))*2.0 +
                                                            hyhz(vLevel) * xixx(vLevel)(iX) * (pressureData(vLevel)(iX + 1, iY, iZ) - pressureData(vLevel)(iX - 1, iY, iZ))*hx(vLevel) +
                                                            hzhx(vLevel) * ety2(vLevel)(iY) * (pressureData(vLevel)(iX, iY + 1, iZ) + pressureData(vLevel)(iX, iY - 1, iZ))*2.0 +
                                                            hzhx(vLevel) * etyy(vLevel)(iY) * (pressureData(vLevel)(iX, iY + 1, iZ) - pressureData(vLevel)(iX, iY - 1, iZ))*hy(vLevel) +
                                                            hxhy(vLevel) * ztz2(vLevel)(iZ) * (pressureData(vLevel)(iX, iY, iZ + 1) + pressureData(vLevel)(iX, iY, iZ - 1))*2.0 +
                                                            hxhy(vLevel) * ztzz(vLevel)(iZ) * (pressureData(vLevel)(iX, iY, iZ + 1) - pressureData(vLevel)(iX, iY, iZ - 1))*hz(vLevel) -
                                                    2.0 * hxhyhz(vLevel) * residualData(vLevel)(iX, iY, iZ))/
                                                   (4.0 * (hyhz(vLevel)*xix2(vLevel)(iX) + hzhx(vLevel)*ety2(vLevel)(iY) + hxhy(vLevel)*ztz2(vLevel)(iZ)));
                    }
                }
            }
        }

        swap(tmpDataArray, pressureData);

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        smothTimeComp += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif
    }

#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
#endif

    imposeBC();

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    smothTimeTran += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif
}

void multigrid_d3::coarsen() {
    real facePoints, edgePoints, vertPoints;

    int i2, j2, k2;
    int pLevel;

    pLevel = vLevel;
    vLevel += 1;

    // Full weighted restriction operation
    for (int i = 0; i <= stagCore(vLevel).ubound(0); i++) {
        i2 = i*2;
        for (int j = 0; j <= stagCore(vLevel).ubound(1); j++) {
            j2 = j*2;
            for (int k = 0; k <= stagCore(vLevel).ubound(2); k++) {
                k2 = k*2;
                facePoints = (tmpDataArray(pLevel)(i2 + 1, j2, k2) + tmpDataArray(pLevel)(i2 - 1, j2, k2) +
                              tmpDataArray(pLevel)(i2, j2 + 1, k2) + tmpDataArray(pLevel)(i2, j2 - 1, k2) +
                              tmpDataArray(pLevel)(i2, j2, k2 + 1) + tmpDataArray(pLevel)(i2, j2, k2 - 1))*0.0625;
                edgePoints = (tmpDataArray(pLevel)(i2 + 1, j2 + 1, k2) + tmpDataArray(pLevel)(i2 + 1, j2 - 1, k2) +
                              tmpDataArray(pLevel)(i2 - 1, j2 - 1, k2) + tmpDataArray(pLevel)(i2 - 1, j2 + 1, k2) +
                              tmpDataArray(pLevel)(i2, j2 + 1, k2 + 1) + tmpDataArray(pLevel)(i2, j2 - 1, k2 + 1) +
                              tmpDataArray(pLevel)(i2, j2 - 1, k2 - 1) + tmpDataArray(pLevel)(i2, j2 + 1, k2 - 1) +
                              tmpDataArray(pLevel)(i2 + 1, j2, k2 + 1) + tmpDataArray(pLevel)(i2 + 1, j2, k2 - 1) +
                              tmpDataArray(pLevel)(i2 - 1, j2, k2 - 1) + tmpDataArray(pLevel)(i2 - 1, j2, k2 + 1))*0.03125;
                vertPoints = (tmpDataArray(pLevel)(i2 + 1, j2 + 1, k2 + 1) +
                              tmpDataArray(pLevel)(i2 + 1, j2 + 1, k2 - 1) +
                              tmpDataArray(pLevel)(i2 + 1, j2 - 1, k2 + 1) +
                              tmpDataArray(pLevel)(i2 - 1, j2 + 1, k2 + 1) +
                              tmpDataArray(pLevel)(i2 + 1, j2 - 1, k2 - 1) +
                              tmpDataArray(pLevel)(i2 - 1, j2 + 1, k2 - 1) +
                              tmpDataArray(pLevel)(i2 - 1, j2 - 1, k2 + 1) +
                              tmpDataArray(pLevel)(i2 - 1, j2 - 1, k2 - 1))*0.015625;

                residualData(vLevel)(i, j, k) = facePoints + edgePoints + vertPoints + tmpDataArray(pLevel)(i2, j2, k2)*0.125;
            }
        }
    }
}


void multigrid_d3::prolong() {
    // Integer values of starting indices, ending indices, and index increments along each direction
    int xSt, xEn, xIn;
    int ySt, yEn, yIn;
    int zSt, zEn, zIn;

    vLevel -= 1;

    // NOTE: Currently interpolating along X first, then Y and finally Z.
    // Test and see if this order is better or the other order, with Z first, then Y and X is better
    // Depending on the order of variables in memory, one of these may give better performance

    // Maybe the below values can be stored in some array instead of recomputing in each prolongation step

    // INTERPOLATE VARIABLE DATA ALONG X-DIRECTION
    xSt = stagCore.lbound(0) + strideValues(vLevel);
    xEn = stagCore.ubound(0) - strideValues(vLevel);
    xIn = strideValues(vLevel+1);

    ySt = stagCore.lbound(1);
    yEn = stagCore.ubound(1);
    yIn = strideValues(vLevel+1);

    zSt = stagCore.lbound(2);
    zEn = stagCore.ubound(2);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iY = ySt; iY <= yEn; iY += yIn) {
            for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
                pressureData(iX, iY, iZ) = (pressureData(iX + strideValues(vLevel), iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/2.0;
            }
        }
    }

    // INTERPOLATE VARIABLE DATA ALONG Y-DIRECTION
    xSt = stagCore.lbound(0);
    xEn = stagCore.ubound(0);
    xIn = strideValues(vLevel);

    ySt = stagCore.lbound(1) + strideValues(vLevel);
    yEn = stagCore.ubound(1) - strideValues(vLevel);
    yIn = strideValues(vLevel+1);

    zSt = stagCore.lbound(2);
    zEn = stagCore.ubound(2);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iY = ySt; iY <= yEn; iY += yIn) {
            for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
                pressureData(iX, iY, iZ) = (pressureData(iX, iY + strideValues(vLevel), iZ) + pressureData(iX, iY - strideValues(vLevel), iZ))/2.0;
            }
        }
    }

    // INTERPOLATE VARIABLE DATA ALONG Z-DIRECTION
    xSt = stagCore.lbound(0);
    xEn = stagCore.ubound(0);
    xIn = strideValues(vLevel);

    ySt = stagCore.lbound(1);
    yEn = stagCore.ubound(1);
    yIn = strideValues(vLevel);

    zSt = stagCore.lbound(2) + strideValues(vLevel);
    zEn = stagCore.ubound(2) - strideValues(vLevel);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iY = ySt; iY <= yEn; iY += yIn) {
            for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
                pressureData(iX, iY, iZ) = (pressureData(iX, iY, iZ + strideValues(vLevel)) + pressureData(iX, iY, iZ - strideValues(vLevel)))/2.0;
            }
        }
    }
}


real multigrid_d3::computeError(const int normOrder) {
    real residualVal = 0.0;

    real tempValue = 0.0;
    real numValLoc = 0.0;
    real denValLoc = 0.0;
    int valCountLoc = 0;

    // This function is called at the finest grid level only.
    // Moreover it called only under the TEST_POISSON flag
    // Hence it is not written to be very fast

    // Problem with Koenig lookup is that when using the function abs with blitz arrays, it automatically computes
    // the absolute of the float values without hitch.
    // When replacing with computing absolute of individual array elements in a loop, ADL chooses a version of
    // abs in the STL which **rounds off** the number.
    // In this case, abs has to be replaced with fabs.
    for (int i = 0; i <= xEnd(0); i += 1) {
        for (int j = 0; j <= yEnd(0); j += 1) {
            for (int k = 0; k <= zEnd(0); k += 1) {
                tempValue = fabs((xix2(0)(i) * (pressureData(0)(i + 1, j, k) - 2.0*pressureData(0)(i, j, k) + pressureData(0)(i - 1, j, k))/(hx(vLevel)*hx(vLevel)) +
                                  xixx(0)(i) * (pressureData(0)(i + 1, j, k) - pressureData(0)(i - 1, j, k))/(2.0*hx(vLevel)) +
                                  ety2(0)(j) * (pressureData(0)(i, j + 1, k) - 2.0*pressureData(0)(i, j, k) + pressureData(0)(i, j - 1, k))/(hy(vLevel)*hy(vLevel)) +
                                  etyy(0)(j) * (pressureData(0)(i, j + 1, k) - pressureData(0)(i, j - 1, k))/(2.0*hy(vLevel)) +
                                  ztz2(0)(k) * (pressureData(0)(i, j, k + 1) - 2.0*pressureData(0)(i, j, k) + pressureData(0)(i, j, k - 1))/(hz(vLevel)*hz(vLevel)) +
                                  ztzz(0)(k) * (pressureData(0)(i, j, k + 1) - pressureData(0)(i, j, k - 1))/(2.0*hz(vLevel))) - residualData(0)(i, j, k));

                switch (normOrder) {
                    case 1:
                        if (tempValue > numValLoc) numValLoc = tempValue;
                        break;
                    case 2:
                        numValLoc += tempValue*tempValue;
                        denValLoc += residualData(0)(i, j, k)*residualData(0)(i, j, k);
                        valCountLoc += 1;
                        break;
                }
            }
        }
    }

    real numValGlo = 0.0;
    real denValGlo = 0.0;
    int valCountGlo = 0;
    switch (normOrder) {
        case 1:
            denValLoc = blitz::max(fabs(residualData(0)));
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
            //residualVal = numValGlo/denValGlo;
            residualVal = numValGlo;
            break;
        case 2:
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&valCountLoc, &valCountGlo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            //residualVal = sqrt(numValGlo/valCountGlo)/sqrt(denValGlo/valCountGlo);
            residualVal = sqrt(numValGlo/valCountGlo);
            break;
    }

    return residualVal;
}


void multigrid_d3::initDirichlet() {
    real xDist, yDist, zDist;

    // Generate the walls as 2D Blitz arrays
    xWall.resize(blitz::TinyVector<int, 2>(stagCore.ubound(1) - stagCore.lbound(1) + 1, stagCore.ubound(2) - stagCore.lbound(2) + 1));
    xWall.reindexSelf(blitz::TinyVector<int, 2>(stagCore.lbound(1), stagCore.lbound(2)));
    xWall = 0.0;

    yWall.resize(blitz::TinyVector<int, 2>(stagCore.ubound(0) - stagCore.lbound(0) + 1, stagCore.ubound(2) - stagCore.lbound(2) + 1));
    yWall.reindexSelf(blitz::TinyVector<int, 2>(stagCore.lbound(0), stagCore.lbound(2)));
    yWall = 0.0;

    zWall.resize(blitz::TinyVector<int, 2>(stagCore.ubound(0) - stagCore.lbound(0) + 1, stagCore.ubound(1) - stagCore.lbound(1) + 1));
    zWall.reindexSelf(blitz::TinyVector<int, 2>(stagCore.lbound(0), stagCore.lbound(1)));
    zWall = 0.0;

    // Compute values at the walls at all V-Cycle depths, using the (r^2)/6 formula
    for (int n=0; n<=inputParams.vcDepth; ++n) {
        // Along X-direction - Left and Right Walls
        xDist = hx(0)*(int(mgSizeArray(localSizeIndex(0))/2) + 1);
        for (int j=stagCore.lbound(1); j<=stagCore.ubound(1); j++) {
            yDist = hy(0)*(j - stagCore.ubound(1)/2);
            for (int k=stagCore.lbound(2); k<=stagCore.ubound(2); k++) {
                zDist = hz(0)*(k - stagCore.ubound(2)/2);

                xWall(j, k) = (xDist*xDist + yDist*yDist + zDist*zDist)/6.0;
            }
        }

        // Along Y-direction - Front and Rear Walls
        yDist = hy(0)*(int(mgSizeArray(localSizeIndex(1))/2) + 1);
        for (int i=stagCore.lbound(0); i<=stagCore.ubound(0); i++) {
            xDist = hx(0)*(i - stagCore.ubound(0)/2);
            for (int k=stagCore.lbound(2); k<=stagCore.ubound(2); k++) {
                zDist = hz(0)*(k - stagCore.ubound(2)/2);

                yWall(i, k) = (xDist*xDist + yDist*yDist + zDist*zDist)/6.0;
            }
        }

        // Along Z-direction - Top and Bottom Walls
        zDist = hz(0)*(int(mgSizeArray(localSizeIndex(2))/2) + 1);
        for (int i=stagCore.lbound(0); i<=stagCore.ubound(0); i++) {
            xDist = hx(0)*(i - stagCore.ubound(0)/2);
            for (int j=stagCore.lbound(1); j<=stagCore.ubound(1); j++) {
                yDist = hy(0)*(j - stagCore.ubound(1)/2);

                zWall(i, j) = (xDist*xDist + yDist*yDist + zDist*zDist)/6.0;
            }
        }
    }
}


void multigrid_d3::imposeBC() {
    if (not inputParams.xPer) {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION ON PRESSURE AT LEFT AND RIGHT WALLS
        if (zeroBC) {
            if (mesh.rankData.xRank == 0) {
                pressureData(vLevel)(-1, all, all) = 0.0;
            }

            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                pressureData(vLevel)(stagCore.ubound(0) + 1, all, all) = 0.0;
            }
        } else {
            if (mesh.rankData.xRank == 0) {
                pressureData(vLevel)(-1, all, all) = xWall(all, all);
            }

            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                pressureData(vLevel)(stagCore.ubound(0) + 1, all, all) = xWall(all, all);
            }
        }
#else
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT LEFT WALL
        if (mesh.rankData.xRank == 0) {
            pressureData(vLevel)(-1, all, all) = pressureData(vLevel)(1, all, all);
        }

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT RIGHT WALL
        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            pressureData(vLevel)(stagCore.ubound(0) + 1, all, all) = pressureData(vLevel)(stagCore.ubound(0) - 1, all, all);
        }
#endif
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (not inputParams.yPer) {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION ON PRESSURE AT FRONT AND BACK WALLS
        if (zeroBC) {
            if (mesh.rankData.yRank == 0) {
                pressureData(vLevel)(all, -1, all) = 0.0;
            }

            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                pressureData(vLevel)(all, stagCore.ubound(1) + 1, all) = 0.0;
            }
        } else {
            if (mesh.rankData.yRank == 0) {
                pressureData(vLevel)(all, -1, all) = yWall(all, all);
            }

            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                pressureData(vLevel)(all, stagCore.ubound(1) + 1, all) = yWall(all, all);
            }
        }
#else
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT FRONT WALL
        if (mesh.rankData.yRank == 0) {
            pressureData(vLevel)(all, -1, all) = pressureData(vLevel)(all, 1, all);
        }

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT BACK WALL
        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
            pressureData(vLevel)(all, stagCore.ubound(1) + 1, all) = pressureData(vLevel)(all, stagCore.ubound(1) - 1, all);
        }
#endif
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (inputParams.zPer) {
        // PERIODIC BOUNDARY CONDITION ON PRESSURE AT BOTTOM WALL
        pressureData(vLevel)(all, all, -1) = pressureData(vLevel)(all, all, stagCore.ubound(2) - 1);

        // PERIODIC BOUNDARY CONDITION ON PRESSURE AT TOP WALL
        pressureData(vLevel)(all, all, stagCore.ubound(2) + 1) = pressureData(vLevel)(all, all, 1);

    } else {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION ON PRESSURE AT BOTTOM AND TOP WALLS
        if (zeroBC) {
            pressureData(vLevel)(all, all, -1) = 0.0;

            pressureData(vLevel)(all, all, stagCore.ubound(2) + 1) = 0.0;
        } else {
            pressureData(vLevel)(all, all, -1) = zWall(all, all);

            pressureData(vLevel)(all, all, stagCore.ubound(2) + 1) = zWall(all, all);
        }
#else
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT BOTTOM WALL
        pressureData(vLevel)(all, all, -1) = pressureData(vLevel)(all, all, 1);

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT TOP WALL
        pressureData(vLevel)(all, all, stagCore.ubound(2) + 1) = pressureData(vLevel)(all, all, stagCore.ubound(2) - 1);
#endif
    }
}
