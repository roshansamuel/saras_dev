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
/*! \file poisson2.cc
 *
 *  \brief Definitions for functions of class poisson for 2D
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
 * \brief   Constructor of the multigrid_d2 class derived from the poisson class
 *
 *          The constructor of the derived multigrid_d2 class frst calls the base poisson class with the arguments passed to it.
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
multigrid_d2::multigrid_d2(const grid &mesh, const parser &solParam): poisson(mesh, solParam) {
    // GET THE localSizeIndex AS IT WILL BE USED TO SET THE FULL AND CORE LIMITS OF THE STAGGERED POINTS
    setLocalSizeIndex();

    // SET THE FULL AND CORE LIMTS SET ABOVE USING THE localSizeIndex VARIBLE SET ABOVE
    setStagBounds();

    // USING THE FULL AND CORE LIMTS SET ABOVE, CREATE ALL Range OBJECTS
    initMeshRanges();

    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // COPY THE STAGGERED GRID DERIVATIVES TO LOCAL ARRAYS
    copyStaggrDerivs();

    // RESIZE AND INITIALIZE NECESSARY DATA-STRUCTURES
    initializeArrays();

    // CREATE THE MPI SUB-ARRAYS NECESSARY TO TRANSFER DATA ACROSS SUB-DOMAINS AT ALL MESH LEVELS
    createMGSubArrays();
}


void multigrid_d2::computeResidual() {
    int iY = 0;

    // Compute Laplacian of the pressure field and subtract it from the RHS of Poisson equation to obtain the residual
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
    for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
        for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
            residualData(iX, iY, iZ) =  inputRHSData(iX, iY, iZ) -
                                       (xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/(hx(vLevel)*hx(vLevel)) +
                                        xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - pressureData(iX - strideValues(vLevel), iY, iZ))/(2.0*hx(vLevel)) +
                                        ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX, iY, iZ - strideValues(vLevel)))/(hz(vLevel)*hz(vLevel)) +
                                        ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - pressureData(iX, iY, iZ - strideValues(vLevel)))/(2.0*hz(vLevel)));
        }
    }
}


void multigrid_d2::smooth(const int smoothCount) {
    iteratorTemp = 0.0;

    for(int n=0; n<smoothCount; n++) {
        // IMPOSE BOUNDARY CONDITION
        imposeBC();

        int iY = 0;
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
            for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
                iteratorTemp(iX, iY, iZ) = (hz2(vLevel) * xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))*2.0 +
                                            hz2(vLevel) * xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - pressureData(iX - strideValues(vLevel), iY, iZ))*hx(vLevel) +
                                            hx2(vLevel) * ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) + pressureData(iX, iY, iZ - strideValues(vLevel)))*2.0 +
                                            hx2(vLevel) * ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - pressureData(iX, iY, iZ - strideValues(vLevel)))*hz(vLevel) -
                                      2.0 * hzhx(vLevel) * residualData(iX, iY, iZ))/
                                    (4.0 * (hz2(vLevel)*xix2(iX) + hx2(vLevel)*ztz2(iZ)));
            }
        }

        swap(iteratorTemp, pressureData);
    }

    imposeBC();
}


void multigrid_d2::coarsen() {
    double facePoints, vertPoints;
    // Integer values of starting indices, ending indices, and index increments along each direction
    int xSt, xEn, xIn;
    int zSt, zEn, zIn;
    int shiftInc;
    int iY = 0;

    imposeBC();

    vLevel += 1;

    xSt = stagCore.lbound(0);
    xEn = stagCore.ubound(0);
    xIn = strideValues(vLevel);

    zSt = stagCore.lbound(2);
    zEn = stagCore.ubound(2);
    zIn = strideValues(vLevel);

    shiftInc = strideValues(vLevel - 1);

    // Full weighted restriction operation
    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
            facePoints = (residualData(iX + shiftInc, iY, iZ) + residualData(iX - shiftInc, iY, iZ) +
                          residualData(iX, iY, iZ + shiftInc) + residualData(iX, iY, iZ - shiftInc))*0.125;
            vertPoints = (residualData(iX + shiftInc, iY, iZ + shiftInc) +
                          residualData(iX + shiftInc, iY, iZ - shiftInc) +
                          residualData(iX - shiftInc, iY, iZ + shiftInc) +
                          residualData(iX - shiftInc, iY, iZ - shiftInc))*0.0625;

            residualData(iX, iY, iZ) = facePoints + vertPoints + residualData(iX, iY, iZ)*0.25;
        }
    }
}


void multigrid_d2::prolong() {
    // Integer values of starting indices, ending indices, and index increments along each direction
    int iY = 0;
    int xSt, xEn, xIn;
    int zSt, zEn, zIn;

    vLevel -= 1;

    // NOTE: Currently interpolating along X first, and then Z.
    // Test and see if this order is better or the other order, with Z first, and then X is better
    // Depending on the order of variables in memory, one of these may give better performance

    // Maybe the below values can be stored in some array instead of recomputing in each prolongation step

    // INTERPOLATE VARIABLE DATA ALONG X-DIRECTION
    xSt = stagCore.lbound(0) + strideValues(vLevel);
    xEn = stagCore.ubound(0) - strideValues(vLevel);
    xIn = strideValues(vLevel+1);

    zSt = stagCore.lbound(2);
    zEn = stagCore.ubound(2);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
            pressureData(iX, iY, iZ) = (pressureData(iX + strideValues(vLevel), iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/2.0;
        }
    }

    // INTERPOLATE VARIABLE DATA ALONG Z-DIRECTION
    xSt = stagCore.lbound(0);
    xEn = stagCore.ubound(0);
    xIn = strideValues(vLevel);

    zSt = stagCore.lbound(2) + strideValues(vLevel);
    zEn = stagCore.ubound(2) - strideValues(vLevel);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
            pressureData(iX, iY, iZ) = (pressureData(iX, iY, iZ + strideValues(vLevel)) + pressureData(iX, iY, iZ - strideValues(vLevel)))/2.0;
        }
    }
}


double multigrid_d2::computeError(const int normOrder) {
    int iY = 0;
    double residualVal = 0.0;

    double tempValue = 0.0;
    double numValLoc = 0.0;
    double denValLoc = 0.0;
    int valCountLoc = 0;

    // Problem with Koenig lookup is that when using the function abs with blitz arrays, it automatically computes
    // the absolute of the float values without hitch.
    // When replacing with computing absolute of individual array elements in a loop, ADL chooses a version of
    // abs in the STL which **rounds off** the number.
    // In this case, abs has to be replaced with fabs.
    for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
        for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
            tempValue = fabs((xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/(hx(vLevel)*hx(vLevel)) +
                              xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - pressureData(iX - strideValues(vLevel), iY, iZ))/(2.0*hx(vLevel)) +
                              ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX, iY, iZ - strideValues(vLevel)))/(hz(vLevel)*hz(vLevel)) +
                              ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - pressureData(iX, iY, iZ - strideValues(vLevel)))/(2.0*hz(vLevel))) - inputRHSData(iX, iY, iZ));

            if (normOrder == 0) {
                if (tempValue > numValLoc) numValLoc = tempValue;
            } else {
                numValLoc += tempValue*tempValue;
                denValLoc += inputRHSData(iX, iY, iZ)*inputRHSData(iX, iY, iZ);
                valCountLoc += 1;
            }
        }
    }

    double numValGlo = 0.0;
    double denValGlo = 0.0;
    int valCountGlo = 0;
    if (normOrder == 0) {
        denValLoc = blitz::max(fabs(inputRHSData));
        MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);
        residualVal = numValGlo/denValGlo;
    } else {
        MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&valCountLoc, &valCountGlo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        residualVal = sqrt(numValGlo/valCountGlo)/sqrt(denValGlo/valCountGlo);
    }
    MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&valCountLoc, &valCountGlo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (mesh.rankData.rank == 0) residualVal = sqrt(numValGlo/valCountGlo)/sqrt(denValGlo/valCountGlo);

    return residualVal;
}


void multigrid_d2::createMGSubArrays() {
    int count, length, stride;

    recvStatus.resize(2);
    recvRequest.resize(2);

    xMGArray.resize(inputParams.vcDepth + 1);
    mgSendLft.resize(inputParams.vcDepth + 1);        mgSendRgt.resize(inputParams.vcDepth + 1);
    mgRecvLft.resize(inputParams.vcDepth + 1);        mgRecvRgt.resize(inputParams.vcDepth + 1);

    for(int i=0; i<=inputParams.vcDepth; i++) {
        // CREATE X_MG_ARRAY DATATYPE
        count = (stagCore.ubound(2) - stagCore.lbound(2))/strideValues(i) + 1;
        length = 1;
        stride = strideValues(i);

        MPI_Type_vector(count, length, stride, MPI_DOUBLE_PRECISION, &xMGArray(i));
        MPI_Type_commit(&xMGArray(i));

        mgSendLft(i) =  strideValues(i), 0, 0;
        mgRecvLft(i) = -strideValues(i), 0, 0;
        mgSendRgt(i) = stagCore.ubound(0) - strideValues(i), 0, 0;
        mgRecvRgt(i) = stagCore.ubound(0) + strideValues(i), 0, 0;

    }
}


void multigrid_d2::imposeBC() {
    updatePads();

    if (not inputParams.xPer) {
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT LEFT WALL
        if (mesh.rankData.xRank == 0) {
            pressureData(-strideValues(vLevel), 0, zMeshRange(vLevel)) = pressureData(strideValues(vLevel), 0, zMeshRange(vLevel));
        }

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT RIGHT WALL
        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            pressureData(stagCore.ubound(0) + strideValues(vLevel), 0, zMeshRange(vLevel)) = pressureData(stagCore.ubound(0) - strideValues(vLevel), 0, zMeshRange(vLevel));
        }
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (inputParams.zPer) {
        // PERIODIC BOUNDARY CONDITION ON PRESSURE AT BOTTOM WALL
        pressureData(xMeshRange(vLevel), 0, -strideValues(vLevel)) = pressureData(xMeshRange(vLevel), 0, stagCore.ubound(2) - strideValues(vLevel));

        // PERIODIC BOUNDARY CONDITION ON PRESSURE AT TOP WALL
        pressureData(xMeshRange(vLevel), 0, stagCore.ubound(2) + strideValues(vLevel)) = pressureData(xMeshRange(vLevel), 0, strideValues(vLevel));

    } else {
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT BOTTOM WALL
        pressureData(xMeshRange(vLevel), 0, -strideValues(vLevel)) = pressureData(xMeshRange(vLevel), 0, strideValues(vLevel));

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT TOP WALL
        pressureData(xMeshRange(vLevel), 0, stagCore.ubound(2) + strideValues(vLevel)) = pressureData(xMeshRange(vLevel), 0, stagCore.ubound(2) - strideValues(vLevel));
    }
}


void multigrid_d2::updatePads() {
    recvRequest = MPI_REQUEST_NULL;

    // TRANSFER DATA FROM NEIGHBOURING CELL TO IMPOSE SUB-DOMAIN BOUNDARY CONDITIONS
    MPI_Irecv(&pressureData(mgRecvLft(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(&pressureData(mgRecvRgt(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));

    MPI_Send(&pressureData(mgSendLft(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(0), 2, MPI_COMM_WORLD);
    MPI_Send(&pressureData(mgSendRgt(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(1), 1, MPI_COMM_WORLD);

    MPI_Waitall(2, recvRequest.dataFirst(), recvStatus.dataFirst());
}


double multigrid_d2::testProlong() {
    int iY = 0;
    vLevel = 0;

    // Fill the residualData array with correct values expected after prolongation
    residualData = 0.0;
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
            residualData(iX, iY, iZ) = (mesh.rankData.rank + 1)*100 + iX*10 + iZ;
        }
    }

    // After going one level down the V-Cycle, populate the pressureData array with values at the corresponding stride
    vLevel += 1;
    pressureData = 0.0;
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
            pressureData(iX, iY, iZ) = (mesh.rankData.rank + 1)*100 + iX*10 + iZ;
        }
    }

    // Perform prolongation
    prolong();

    pressureData -= residualData;

    return blitz::max(fabs(pressureData));
}


double multigrid_d2::testTransfer() {
    double maxVal = 0.0;

    int iY = 0;
    vLevel = 0;

    pressureData = 0.0;
    residualData = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += 1) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += 1) {
            pressureData(iX, iY, iZ) = (mesh.rankData.rank + 1)*100 + iX*10 + iZ;
            residualData(iX, iY, iZ) = pressureData(iX, iY, iZ);
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int iX = 0; iX <= inputParams.vcDepth; iX++) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(iX)) {
            residualData(-strideValues(iX), iY, iZ) = (mesh.rankData.nearRanks(0) + 1)*100 + (stagCore.ubound(0) - strideValues(iX))*10 + iZ;
            residualData(stagCore.ubound(0) + strideValues(iX), iY, iZ) = (mesh.rankData.nearRanks(1) + 1)*100 + strideValues(iX)*10 + iZ;
        }
    }

    for (int i=0; i<=inputParams.vcDepth; i++) {
        updatePads();
        vLevel += 1;
    }

    pressureData -= residualData;

    for (int iX = pressureData.lbound(0); iX <= pressureData.ubound(0); iX += 1) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += 1) {
            if (abs(pressureData(iX, iY, iZ)) > maxVal) {
                maxVal = abs(pressureData(iX, iY, iZ));
            }
        }
    }

    return maxVal;
}


double multigrid_d2::testPeriodic() {
    int iY = 0;
    double xCoord = 0.0;
    double zCoord = 0.0;

    vLevel = 0;

    pressureData = 0.0;
    residualData = 0.0;

    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
            pressureData(iX, iY, iZ) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                       cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            residualData(iX, iY, iZ) = pressureData(iX, iY, iZ);
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int iX = 0; iX <= inputParams.vcDepth; iX++) {
        for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(iX)) {
            xCoord = mesh.xStaggr(stagCore.lbound(0)) - (mesh.xStaggr(stagCore.lbound(0) + strideValues(iX)) - mesh.xStaggr(stagCore.lbound(0)));
            residualData(stagCore.lbound(0) - strideValues(iX), iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                                          cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

            xCoord = mesh.xStaggr(stagCore.ubound(0)) + (mesh.xStaggr(stagCore.ubound(0)) - mesh.xStaggr(stagCore.ubound(0) - strideValues(iX)));
            residualData(stagCore.ubound(0) + strideValues(iX), iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                                          cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
        }
    }

    for (int iZ = 0; iZ <= inputParams.vcDepth; iZ++) {
        for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(iZ)) {
            zCoord = mesh.zStaggr(stagCore.lbound(2)) - (mesh.zStaggr(stagCore.lbound(2) + strideValues(iZ)) - mesh.zStaggr(stagCore.lbound(2)));
            residualData(iX, iY, stagCore.lbound(2) - strideValues(iZ)) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                          cos(2.0*M_PI*zCoord/mesh.zLen);

            zCoord = mesh.zStaggr(stagCore.ubound(2)) + (mesh.zStaggr(stagCore.ubound(2)) - mesh.zStaggr(stagCore.ubound(2) - strideValues(iZ)));
            residualData(iX, iY, stagCore.ubound(2) + strideValues(iZ)) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                          cos(2.0*M_PI*zCoord/mesh.zLen);
        }
    }

    for (int i=0; i<=inputParams.vcDepth; i++) {
        imposeBC();
        vLevel += 1;
    }

    pressureData -= residualData;

    return blitz::max(fabs(pressureData));
}
