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
 *          Most of the array initializations are performed in the constructor of the base class.
 *          Initializations unique to 3D implementation of multi-grid solver are then done here.
 *          Primarily this involves generating the MPI datatypes for data transfer between sub-domain boundaries.
 *          When testing the Poisson solver with Dirichlet BC, the analytical 3D solution for the test is also calculated here
 *          through a call to initDirichlet() function.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
multigrid_d3::multigrid_d3(const grid &mesh, const parser &solParam): poisson(mesh, solParam) {
    // CREATE THE MPI SUB-ARRAYS NECESSARY TO TRANSFER DATA ACROSS SUB-DOMAINS AT ALL MESH LEVELS
    createMGSubArrays();

    // INITIALIZE DIRICHLET BCs WHEN TESTING THE POISSON SOLVER
#ifdef TEST_POISSON
    initDirichlet();
#endif
}


void multigrid_d3::computeResidual() {
    tmpDataArray(vLevel) = 0.0;

    // Compute Laplacian of the pressure field and subtract it from the RHS of Poisson equation to obtain the residual
    // This residual is temporarily stored into tmpDataArray, from which it will be coarsened into residualData array.
    // Needed update: Substitute the below OpenMP parallel loop with vectorized Blitz operation and check for speed increase.
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int j = 0; j <= yEnd(vLevel); ++j) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                tmpDataArray(vLevel)(i, j, k) =  residualData(vLevel)(i, j, k) -
                             (xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - 2.0*pressureData(vLevel)(i, j, k) + pressureData(vLevel)(i - 1, j, k))/(hx(vLevel)*hx(vLevel)) +
                              xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - pressureData(vLevel)(i - 1, j, k))/(2.0*hx(vLevel)) +
                              ety2(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - 2.0*pressureData(vLevel)(i, j, k) + pressureData(vLevel)(i, j - 1, k))/(hy(vLevel)*hy(vLevel)) +
                              etyy(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - pressureData(vLevel)(i, j - 1, k))/(2.0*hy(vLevel)) +
                              ztz2(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - 2.0*pressureData(vLevel)(i, j, k) + pressureData(vLevel)(i, j, k - 1))/(hz(vLevel)*hz(vLevel)) +
                              ztzz(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - pressureData(vLevel)(i, j, k - 1))/(2.0*hz(vLevel)));
            }
        }
    }

    updatePads(tmpDataArray);
}


void multigrid_d3::smooth(const int smoothCount) {
    tmpDataArray(vLevel) = 0.0;

    for(int n=0; n<smoothCount; ++n) {
        imposeBC();

        // WARNING: When using the gauss-seidel smoothing as written below, the edges of interior sub-domains after MPI decomposition will not have the updated values
        // As a result, the serial and parallel results will not match when using gauss-seidel smoothing
        if (inputParams.gsSmooth) {
            // GAUSS-SEIDEL ITERATIVE SMOOTHING
            for (int i = 0; i <= xEnd(vLevel); ++i) {
                for (int j = 0; j <= yEnd(vLevel); ++j) {
                    for (int k = 0; k <= zEnd(vLevel); ++k) {
                        pressureData(vLevel)(i, j, k) = (hyhz(vLevel) * xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) + pressureData(vLevel)(i - 1, j, k))*2.0 +
                                                         hyhz(vLevel) * xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - pressureData(vLevel)(i - 1, j, k))*hx(vLevel) +
                                                         hzhx(vLevel) * ety2(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) + pressureData(vLevel)(i, j - 1, k))*2.0 +
                                                         hzhx(vLevel) * etyy(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - pressureData(vLevel)(i, j - 1, k))*hy(vLevel) +
                                                         hxhy(vLevel) * ztz2(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) + pressureData(vLevel)(i, j, k - 1))*2.0 +
                                                         hxhy(vLevel) * ztzz(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - pressureData(vLevel)(i, j, k - 1))*hz(vLevel) -
                                                 2.0 * hxhyhz(vLevel) * residualData(vLevel)(i, j, k))/
                                                (4.0 * (hyhz(vLevel)*xix2(vLevel)(i) + hzhx(vLevel)*ety2(vLevel)(j) + hxhy(vLevel)*ztz2(vLevel)(k)));
                    }
                }
            }
        } else {
            // JACOBI ITERATIVE SMOOTHING
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
            for (int i = 0; i <= xEnd(vLevel); ++i) {
                for (int j = 0; j <= yEnd(vLevel); ++j) {
                    for (int k = 0; k <= zEnd(vLevel); ++k) {
                        tmpDataArray(vLevel)(i, j, k) = (hyhz(vLevel) * xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) + pressureData(vLevel)(i - 1, j, k))*2.0 +
                                                         hyhz(vLevel) * xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - pressureData(vLevel)(i - 1, j, k))*hx(vLevel) +
                                                         hzhx(vLevel) * ety2(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) + pressureData(vLevel)(i, j - 1, k))*2.0 +
                                                         hzhx(vLevel) * etyy(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - pressureData(vLevel)(i, j - 1, k))*hy(vLevel) +
                                                         hxhy(vLevel) * ztz2(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) + pressureData(vLevel)(i, j, k - 1))*2.0 +
                                                         hxhy(vLevel) * ztzz(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - pressureData(vLevel)(i, j, k - 1))*hz(vLevel) -
                                                 2.0 * hxhyhz(vLevel) * residualData(vLevel)(i, j, k))/
                                                (4.0 * (hyhz(vLevel)*xix2(vLevel)(i) + hzhx(vLevel)*ety2(vLevel)(j) + hxhy(vLevel)*ztz2(vLevel)(k)));
                    }
                }
            }

            swap(tmpDataArray, pressureData);
        }
    }

    imposeBC();
}


void multigrid_d3::solve() {
    int iterCount = 0;
    real tempValue, localMax, globalMax;

    while (true) {
        imposeBC();

        // GAUSS-SEIDEL ITERATIVE SOLVER
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int j = 0; j <= yEnd(vLevel); ++j) {
                for (int k = 0; k <= zEnd(vLevel); ++k) {
                    pressureData(vLevel)(i, j, k) = (hyhz(vLevel) * xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) + pressureData(vLevel)(i - 1, j, k))*2.0 +
                                                     hyhz(vLevel) * xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - pressureData(vLevel)(i - 1, j, k))*hx(vLevel) +
                                                     hzhx(vLevel) * ety2(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) + pressureData(vLevel)(i, j - 1, k))*2.0 +
                                                     hzhx(vLevel) * etyy(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - pressureData(vLevel)(i, j - 1, k))*hy(vLevel) +
                                                     hxhy(vLevel) * ztz2(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) + pressureData(vLevel)(i, j, k - 1))*2.0 +
                                                     hxhy(vLevel) * ztzz(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - pressureData(vLevel)(i, j, k - 1))*hz(vLevel) -
                                             2.0 * hxhyhz(vLevel) * residualData(vLevel)(i, j, k))/
                                             (4.0 * (hyhz(vLevel) * xix2(vLevel)(i) + hzhx(vLevel) * ety2(vLevel)(j) + hxhy(vLevel) * ztz2(vLevel)(k)));
                }
            }
        }

        tempValue = 0.0;
        localMax = -1.0e-10;
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int j = 0; j <= yEnd(vLevel); ++j) {
                for (int k = 0; k <= zEnd(vLevel); ++k) {
                    tempValue = fabs(residualData(vLevel)(i, j, k) -
                               (xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - 2.0*pressureData(vLevel)(i, j, k) + pressureData(vLevel)(i - 1, j, k))/(hx(vLevel)*hx(vLevel)) +
                                xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, j, k) - pressureData(vLevel)(i - 1, j, k))/(2.0*hx(vLevel)) +
                                ety2(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - 2.0*pressureData(vLevel)(i, j, k) + pressureData(vLevel)(i, j - 1, k))/(hy(vLevel)*hy(vLevel)) +
                                etyy(vLevel)(j) * (pressureData(vLevel)(i, j + 1, k) - pressureData(vLevel)(i, j - 1, k))/(2.0*hy(vLevel)) +
                                ztz2(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - 2.0*pressureData(vLevel)(i, j, k) + pressureData(vLevel)(i, j, k - 1))/(hz(vLevel)*hz(vLevel)) +
                                ztzz(vLevel)(k) * (pressureData(vLevel)(i, j, k + 1) - pressureData(vLevel)(i, j, k - 1))/(2.0*hz(vLevel))));

                    if (tempValue > localMax) {
                        localMax = tempValue;
                    }
                }
            }
        }

        MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

        if (globalMax < inputParams.mgTolerance) {
            break;
        }

        iterCount += 1;
        if (iterCount > maxCount) {
            if (inputParams.printResidual) if (mesh.rankData.rank == 0) std::cout << "WARNING: Iterations for solution at coarsest level not converging." << std::endl;
            break;
        }
    }

    imposeBC();
}


void multigrid_d3::coarsen() {
    real facePoints, edgePoints, vertPoints, centrPoint;
    real v000, v001, v010, v100, v011, v101, v110, v111, vTot;

    int i2, j2, k2;
    int pLevel;

    pLevel = vLevel;
    vLevel += 1;

    // Full weighted restriction operation
    // The residual computed at previous vLevel is stored in tmpDataArray.
    // This data is read for coarsening and written into residualData array.

    /*
     * According to An Introduction to Multigrid Methods by P. Wesseling, Page 64 (Sec 5.2),
     * Restriction can be performed at the edges and corners using the same stencil as in the bulk,
     * But by assuming that values of the field outside the domain are all 0.
     * Hence no special treatments at the corners and edges are needed.
     */

    for (int i = 0; i <= xEnd(vLevel); ++i) {
        i2 = i*2;
        for (int j = 0; j <= yEnd(vLevel); ++j) {
            j2 = j*2;
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                k2 = k*2;

                // Volumes of the 8 cubiodal sub-sections into which the coarse grid point divides the fine grid cell.
                // These values reduce to fractions of powers of 2 for uniform grids.
                // The total volume of the coarse grid cell is also calculated.
                v000 = (nuX(vLevel)(i) - nuX(pLevel)(i2 - 1))*(nuY(vLevel)(j) - nuY(pLevel)(j2 - 1))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2 - 1));
                v100 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(vLevel)(j) - nuY(pLevel)(j2 - 1))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2 - 1));
                v010 = (nuX(vLevel)(i) - nuX(pLevel)(i2 - 1))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2 - 1));
                v001 = (nuX(vLevel)(i) - nuX(pLevel)(i2 - 1))*(nuY(vLevel)(j) - nuY(pLevel)(j2 - 1))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                v011 = (nuX(vLevel)(i) - nuX(pLevel)(i2 - 1))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                v101 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(vLevel)(j) - nuY(pLevel)(j2 - 1))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                v110 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2 - 1));
                v111 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                vTot = (nuX(pLevel)(i2 + 1) - nuX(pLevel)(i2 - 1))*(nuY(pLevel)(j2 + 1) - nuY(pLevel)(j2 - 1))*(nuZ(pLevel)(k2 + 1) - nuZ(pLevel)(k2 - 1));

                // Full-weighted restriction for non-uniform grids
                facePoints = ((v000 + v001 + v010 + v011)*tmpDataArray(pLevel)(i2 - 1, j2, k2) + (v100 + v101 + v110 + v111)*tmpDataArray(pLevel)(i2 + 1, j2, k2) +
                              (v000 + v001 + v100 + v101)*tmpDataArray(pLevel)(i2, j2 - 1, k2) + (v010 + v011 + v110 + v111)*tmpDataArray(pLevel)(i2, j2 + 1, k2) +
                              (v000 + v100 + v010 + v110)*tmpDataArray(pLevel)(i2, j2, k2 - 1) + (v001 + v101 + v011 + v111)*tmpDataArray(pLevel)(i2, j2, k2 + 1));
                edgePoints = ((v110 + v111)*tmpDataArray(pLevel)(i2 + 1, j2 + 1, k2) + (v100 + v101)*tmpDataArray(pLevel)(i2 + 1, j2 - 1, k2) +
                              (v000 + v001)*tmpDataArray(pLevel)(i2 - 1, j2 - 1, k2) + (v010 + v011)*tmpDataArray(pLevel)(i2 - 1, j2 + 1, k2) +
                              (v011 + v111)*tmpDataArray(pLevel)(i2, j2 + 1, k2 + 1) + (v001 + v101)*tmpDataArray(pLevel)(i2, j2 - 1, k2 + 1) +
                              (v000 + v100)*tmpDataArray(pLevel)(i2, j2 - 1, k2 - 1) + (v010 + v110)*tmpDataArray(pLevel)(i2, j2 + 1, k2 - 1) +
                              (v101 + v111)*tmpDataArray(pLevel)(i2 + 1, j2, k2 + 1) + (v100 + v110)*tmpDataArray(pLevel)(i2 + 1, j2, k2 - 1) +
                              (v000 + v010)*tmpDataArray(pLevel)(i2 - 1, j2, k2 - 1) + (v001 + v011)*tmpDataArray(pLevel)(i2 - 1, j2, k2 + 1));
                vertPoints = (v111*tmpDataArray(pLevel)(i2 + 1, j2 + 1, k2 + 1) +
                              v110*tmpDataArray(pLevel)(i2 + 1, j2 + 1, k2 - 1) +
                              v101*tmpDataArray(pLevel)(i2 + 1, j2 - 1, k2 + 1) +
                              v011*tmpDataArray(pLevel)(i2 - 1, j2 + 1, k2 + 1) +
                              v100*tmpDataArray(pLevel)(i2 + 1, j2 - 1, k2 - 1) +
                              v010*tmpDataArray(pLevel)(i2 - 1, j2 + 1, k2 - 1) +
                              v001*tmpDataArray(pLevel)(i2 - 1, j2 - 1, k2 + 1) +
                              v000*tmpDataArray(pLevel)(i2 - 1, j2 - 1, k2 - 1));
                centrPoint = (vTot*tmpDataArray(pLevel)(i2, j2, k2));

                residualData(vLevel)(i, j, k) = (facePoints + edgePoints + vertPoints + centrPoint)/(vTot*8.0);
            }
        }
    }
}


void multigrid_d3::prolong() {
    int pLevel;
    int i2, j2, k2;

    pLevel = vLevel;
    vLevel -= 1;

    pressureData(vLevel) = 0.0;

    // This method of using loops with conditional expressions was much faster (nearly 1.5 times)
    // than the old method of using 3 vectorized sweeps with Blitz Range objects on uniform grids
    // without additional computational work of calculating interpolation coefficients on
    // non-uniform grids.
    real c0, c1, c2, c3, c4, c5, c6, c7, c8;
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        i2 = i/2;
        if (isOdd(i)) {
            for (int j = 0; j <= yEnd(vLevel); ++j) {
                j2 = j/2;
                if (isOdd(j)) {
                    for (int k = 0; k <= zEnd(vLevel); ++k) {
                        k2 = k/2;
                        if (isOdd(k)) {     // i j and k are odd
                            c0 = (nuX(pLevel)(i2 + 1) - nuX(pLevel)(i2))*(nuY(pLevel)(j2 + 1) - nuY(pLevel)(j2))*(nuZ(pLevel)(k2 + 1) - nuZ(pLevel)(k2));
                            c1 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c2 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c3 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(vLevel)(j) - nuY(pLevel)(j2))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c4 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));
                            c5 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuY(vLevel)(j) - nuY(pLevel)(j2))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c6 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));
                            c7 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(vLevel)(j) - nuY(pLevel)(j2))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));
                            c8 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuY(vLevel)(j) - nuY(pLevel)(j2))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) +
                                                             c2*pressureData(pLevel)(i2 + 1, j2, k2) + c3*pressureData(pLevel)(i2, j2 + 1, k2) + c4*pressureData(pLevel)(i2, j2, k2 + 1) +
                                                             c5*pressureData(pLevel)(i2 + 1, j2 + 1, k2) + c6*pressureData(pLevel)(i2 + 1, j2, k2 + 1) + c7*pressureData(pLevel)(i2, j2 + 1, k2 + 1) +
                                                             c8*pressureData(pLevel)(i2 + 1, j2 + 1, k2 + 1))/c0;

                        } else {            // i and j are odd, but k is even
                            c0 = (nuX(pLevel)(i2 + 1) - nuX(pLevel)(i2))*(nuY(pLevel)(j2 + 1) - nuY(pLevel)(j2));
                            c1 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j));
                            c2 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuY(vLevel)(j) - nuY(pLevel)(j2));
                            c3 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuY(pLevel)(j2 + 1) - nuY(vLevel)(j));
                            c4 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuY(vLevel)(j) - nuY(pLevel)(j2));

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) + c2*pressureData(pLevel)(i2, j2 + 1, k2) +
                                                             c3*pressureData(pLevel)(i2 + 1, j2, k2) + c4*pressureData(pLevel)(i2 + 1, j2 + 1, k2))/c0;
                        }
                    }
                } else {
                    for (int k = 0; k <= zEnd(vLevel); ++k) {
                        k2 = k/2;
                        if (isOdd(k)) {     // i and k are odd, but j is even
                            c0 = (nuX(pLevel)(i2 + 1) - nuX(pLevel)(i2))*(nuZ(pLevel)(k2 + 1) - nuZ(pLevel)(k2));
                            c1 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c2 = (nuX(pLevel)(i2 + 1) - nuX(vLevel)(i))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));
                            c3 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c4 = (nuX(vLevel)(i) - nuX(pLevel)(i2))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) + c2*pressureData(pLevel)(i2, j2, k2 + 1) +
                                                             c3*pressureData(pLevel)(i2 + 1, j2, k2) + c4*pressureData(pLevel)(i2 + 1, j2, k2 + 1))/c0;

                        } else {            // i is odd, but j and k are even
                            c0 = nuX(pLevel)(i2 + 1) - nuX(pLevel)(i2);
                            c1 = nuX(pLevel)(i2 + 1) - nuX(vLevel)(i);
                            c2 = nuX(vLevel)(i) - nuX(pLevel)(i2);

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) + c2*pressureData(pLevel)(i2 + 1, j2, k2))/c0;
                        }
                    }
                }
            }
        } else {
            for (int j = 0; j <= yEnd(vLevel); j++) {
                j2 = j/2;
                if (isOdd(j)) {
                    for (int k = 0; k <= zEnd(vLevel); ++k) {
                        k2 = k/2;
                        if (isOdd(k)) {     // i is even, but j and k are odd
                            c0 = (nuY(pLevel)(j2 + 1) - nuY(pLevel)(j2))*(nuZ(pLevel)(k2 + 1) - nuZ(pLevel)(k2));
                            c1 = (nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c2 = (nuY(pLevel)(j2 + 1) - nuY(vLevel)(j))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));
                            c3 = (nuY(vLevel)(j) - nuY(pLevel)(j2))*(nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k));
                            c4 = (nuY(vLevel)(j) - nuY(pLevel)(j2))*(nuZ(vLevel)(k) - nuZ(pLevel)(k2));

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) + c2*pressureData(pLevel)(i2, j2, k2 + 1) +
                                                             c3*pressureData(pLevel)(i2, j2 + 1, k2) + c4*pressureData(pLevel)(i2, j2 + 1, k2 + 1))/c0;

                        } else {            // i and k are even, but j is odd
                            c0 = nuY(pLevel)(j2 + 1) - nuY(pLevel)(j2);
                            c1 = nuY(pLevel)(j2 + 1) - nuY(vLevel)(j);
                            c2 = nuY(vLevel)(j) - nuY(pLevel)(j2);

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) + c2*pressureData(pLevel)(i2, j2 + 1, k2))/c0;
                        }
                    }
                } else {
                    for (int k = 0; k <= zEnd(vLevel); ++k) {
                        k2 = k/2;
                        if (isOdd(k)) {     // i and j are even, but k is odd
                            c0 = nuZ(pLevel)(k2 + 1) - nuZ(pLevel)(k2);
                            c1 = nuZ(pLevel)(k2 + 1) - nuZ(vLevel)(k);
                            c2 = nuZ(vLevel)(k) - nuZ(pLevel)(k2);

                            pressureData(vLevel)(i, j, k) = (c1*pressureData(pLevel)(i2, j2, k2) + c2*pressureData(pLevel)(i2, j2, k2 + 1))/c0;

                        } else {            // i j and k are even
                            pressureData(vLevel)(i, j, k) = pressureData(pLevel)(i2, j2, k2);
                        }
                    }
                }
            }
        }
    }
}


real multigrid_d3::computeError(const int normOrder) {
    real residualVal = 0.0;
    real numValLoc = 0.0;
    real denValLoc = 0.0;
    real tempNum = 0.0;
    real tempDen = 0.0;

    // This function is called at the finest grid level only.

    // Problem with Koenig lookup is that when using the function abs with blitz arrays, it automatically computes
    // the absolute of the float values without hitch.
    // When replacing with computing absolute of individual array elements in a loop, ADL chooses a version of
    // abs in the STL which **rounds off** the number.
    // In this case, abs has to be replaced with fabs.
    for (int i = 0; i <= xEnd(0); ++i) {
        for (int j = 0; j <= yEnd(0); ++j) {
            for (int k = 0; k <= zEnd(0); ++k) {
                tempNum = fabs((xix2(0)(i) * (pressureData(0)(i + 1, j, k) - 2.0*pressureData(0)(i, j, k) + pressureData(0)(i - 1, j, k))/(hx(vLevel)*hx(vLevel)) +
                                xixx(0)(i) * (pressureData(0)(i + 1, j, k) - pressureData(0)(i - 1, j, k))/(2.0*hx(vLevel)) +
                                ety2(0)(j) * (pressureData(0)(i, j + 1, k) - 2.0*pressureData(0)(i, j, k) + pressureData(0)(i, j - 1, k))/(hy(vLevel)*hy(vLevel)) +
                                etyy(0)(j) * (pressureData(0)(i, j + 1, k) - pressureData(0)(i, j - 1, k))/(2.0*hy(vLevel)) +
                                ztz2(0)(k) * (pressureData(0)(i, j, k + 1) - 2.0*pressureData(0)(i, j, k) + pressureData(0)(i, j, k - 1))/(hz(vLevel)*hz(vLevel)) +
                                ztzz(0)(k) * (pressureData(0)(i, j, k + 1) - pressureData(0)(i, j, k - 1))/(2.0*hz(vLevel))) - residualData(0)(i, j, k));

                tempDen = fabs(residualData(0)(i, j, k));

                switch (normOrder) {
                    case 0: // L-Infinity Norm
                        if (tempNum > numValLoc) numValLoc = tempNum;
                        if (tempDen > denValLoc) denValLoc = tempDen;
                        break;
                    case 1: // L-1 Norm
                        numValLoc += tempNum;
                        denValLoc += tempDen;
                        break;
                    case 2: // L-2 Norm
                        numValLoc += tempNum*tempNum;
                        denValLoc += tempDen*tempDen;
                        break;
                }
            }
        }
    }

    real numValGlo = 0.0;
    real denValGlo = 0.0;
    switch (normOrder) {
        case 0:     // L-Infinity Norm
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

            if (denValGlo) {
                residualVal = numValGlo/denValGlo;
            } else {
                residualVal = numValGlo;
            }
            break;
        case 1:     // L-1 Norm
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);

            if (denValGlo) {
                residualVal = numValGlo/denValGlo;
            } else {
                residualVal = numValGlo/pointCount;
            }
            break;
        case 2:     // L-2 Norm
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);

            if (denValGlo) {
                residualVal = sqrt(numValGlo/pointCount)/sqrt(denValGlo/pointCount);
            } else {
                residualVal = sqrt(numValGlo/pointCount);
            }
            break;
    }

    return residualVal;
}


void multigrid_d3::createMGSubArrays() {
    int count, length, stride;

    recvStatus.resize(4);
    recvRequest.resize(4);

    xMGArray.resize(inputParams.vcDepth + 1);
    yMGArray.resize(inputParams.vcDepth + 1);
    zMGArray.resize(inputParams.vcDepth + 1);

    mgSendLft.resize(inputParams.vcDepth + 1);          mgSendRgt.resize(inputParams.vcDepth + 1);
    mgRecvLft.resize(inputParams.vcDepth + 1);          mgRecvRgt.resize(inputParams.vcDepth + 1);
    mgSendFrn.resize(inputParams.vcDepth + 1);          mgSendBak.resize(inputParams.vcDepth + 1);
    mgRecvFrn.resize(inputParams.vcDepth + 1);          mgRecvBak.resize(inputParams.vcDepth + 1);

    mgSendLftFrn.resize(inputParams.vcDepth + 1);       mgSendRgtBak.resize(inputParams.vcDepth + 1);
    mgRecvLftFrn.resize(inputParams.vcDepth + 1);       mgRecvRgtBak.resize(inputParams.vcDepth + 1);
    mgSendRgtFrn.resize(inputParams.vcDepth + 1);       mgSendLftBak.resize(inputParams.vcDepth + 1);
    mgRecvRgtFrn.resize(inputParams.vcDepth + 1);       mgRecvLftBak.resize(inputParams.vcDepth + 1);

    /***************************************************************************************************
    * Previously xMGArray and yMGArray were defined only if npX > 1 or npY > 1 respectively.
    * This condition remained as a hidden bug in the code for the long time
    * Because for periodic cases, it was implicitly assumed that periodic data transfer will serve
    * But for a sequential case with npX = 1 and npY = 1, this transfer will not happen
    * Now xMGArray and yMGArray are defined irrespective of npX and npY
    \**************************************************************************************************/

    for(int n=0; n<=inputParams.vcDepth; n++) {
        // CREATE X_MG_ARRAY DATATYPE
        count = (stagFull(n).ubound(2) + 2)*(stagFull(n).ubound(1) + 2);

        MPI_Type_contiguous(count, MPI_FP_REAL, &xMGArray(n));
        MPI_Type_commit(&xMGArray(n));

        // CREATE Y_MG_ARRAY DATATYPE
        count = stagFull(n).ubound(0) + 2;
        length = stagFull(n).ubound(2) + 2;
        stride = (stagFull(n).ubound(2) + 2)*(stagFull(n).ubound(1) + 2);

        MPI_Type_vector(count, length, stride, MPI_FP_REAL, &yMGArray(n));
        MPI_Type_commit(&yMGArray(n));

        // CREATE Z_MG_ARRAY DATATYPE - FOR DATA-TRANSFER ACROSS SUB-DOMAIN EDGES
        count = stagFull(n).ubound(2) + 2;

        MPI_Type_contiguous(count, MPI_FP_REAL, &zMGArray(n));
        MPI_Type_commit(&zMGArray(n));

        // SET STARTING INDICES OF MEMORY LOCATIONS FROM WHERE TO READ (SEND) AND WRITE (RECEIVE) DATA
        mgSendLft(n) =  1, -1, -1;
        mgRecvLft(n) = -1, -1, -1;
        mgSendRgt(n) = stagCore(n).ubound(0) - 1, -1, -1;
        mgRecvRgt(n) = stagCore(n).ubound(0) + 1, -1, -1;

        mgSendFrn(n) = -1,  1, -1;
        mgRecvFrn(n) = -1, -1, -1;
        mgSendBak(n) = -1, stagCore(n).ubound(1) - 1, -1;
        mgRecvBak(n) = -1, stagCore(n).ubound(1) + 1, -1;

        mgSendLftFrn(n) =  1,  1, -1;
        mgRecvLftFrn(n) = -1, -1, -1;
        mgSendRgtBak(n) = stagCore(n).ubound(0) - 1, stagCore(n).ubound(1) - 1, -1;
        mgRecvRgtBak(n) = stagCore(n).ubound(0) + 1, stagCore(n).ubound(1) + 1, -1;

        mgSendRgtFrn(n) = stagCore(n).ubound(0) - 1,  1, -1;
        mgRecvRgtFrn(n) = stagCore(n).ubound(0) + 1, -1, -1;
        mgSendLftBak(n) =  1, stagCore(n).ubound(1) - 1, -1;
        mgRecvLftBak(n) = -1, stagCore(n).ubound(1) + 1, -1;
    }
}


void multigrid_d3::initDirichlet() {
    real xDist, yDist, zDist;

    // Generate the walls as 2D Blitz arrays
    xWall.resize(stagFull(0).ubound(1) - stagFull(0).lbound(1) + 1, stagFull(0).ubound(2) - stagFull(0).lbound(2) + 1);
    xWall.reindexSelf(blitz::TinyVector<int, 2>(stagFull(0).lbound(1), stagFull(0).lbound(2)));
    xWall = 0.0;

    yWall.resize(stagFull(0).ubound(0) - stagFull(0).lbound(0) + 1, stagFull(0).ubound(2) - stagFull(0).lbound(2) + 1);
    yWall.reindexSelf(blitz::TinyVector<int, 2>(stagFull(0).lbound(0), stagFull(0).lbound(2)));
    yWall = 0.0;

    zWall.resize(stagFull(0).ubound(0) - stagFull(0).lbound(0) + 1, stagFull(0).ubound(1) - stagFull(0).lbound(1) + 1);
    zWall.reindexSelf(blitz::TinyVector<int, 2>(stagFull(0).lbound(0), stagFull(0).lbound(1)));
    zWall = 0.0;

    // Compute values at the walls using the (r^2)/6 formula
    // Along X-direction - Left and Right Walls
    xDist = mesh.inputParams.Lx/2.0;
    for (int j=0; j<=stagCore(0).ubound(1); ++j) {
        yDist = mesh.y(j) - 0.5;
        for (int k=0; k<=stagCore(0).ubound(2); ++k) {
            zDist = mesh.z(k) - 0.5;

            xWall(j, k) = (xDist*xDist + yDist*yDist + zDist*zDist)/6.0;
        }
    }

    // Along Y-direction - Front and Rear Walls
    yDist = mesh.inputParams.Ly/2.0;
    for (int i=0; i<=stagCore(0).ubound(0); ++i) {
        xDist = mesh.x(i) - 0.5;
        for (int k=0; k<=stagCore(0).ubound(2); ++k) {
            zDist = mesh.z(k) - 0.5;

            yWall(i, k) = (xDist*xDist + yDist*yDist + zDist*zDist)/6.0;
        }
    }

    // Along Z-direction - Top and Bottom Walls
    zDist = mesh.inputParams.Lz/2.0;
    for (int i=0; i<=stagCore(0).ubound(0); ++i) {
        xDist = mesh.x(i) - 0.5;
        for (int j=0; j<=stagCore(0).ubound(1); ++j) {
            yDist = mesh.y(j) - 0.5;

            zWall(i, j) = (xDist*xDist + yDist*yDist + zDist*zDist)/6.0;
        }
    }
}


void multigrid_d3::imposeBC() {
    // FOR PARALLEL RUNS, FIRST UPDATE GHOST POINTS OF MPI SUB-DOMAINS
    updatePads(pressureData);

    if (not inputParams.xPer) {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION AT LEFT AND RIGHT WALLS
        if (zeroBC) {
            if (mesh.rankData.xRank == 0) {
                pressureData(vLevel)(-1, all, all) = -pressureData(vLevel)(1, all, all);
            }

            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, all, all) = -pressureData(vLevel)(stagCore(vLevel).ubound(0) - 1, all, all);
            }
        } else {
            if (mesh.rankData.xRank == 0) {
                pressureData(vLevel)(-1, all, all) = 2.0*xWall(all, all) - pressureData(vLevel)(1, all, all);
            }

            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, all, all) = 2.0*xWall(all, all) - pressureData(vLevel)(stagCore(vLevel).ubound(0) - 1, all, all);
            }
        }
#else
        // NEUMANN BOUNDARY CONDITION AT LEFT AND RIGHT WALLS
        if (mesh.rankData.xRank == 0) {
            pressureData(vLevel)(-1, all, all) = pressureData(vLevel)(1, all, all);
        }

        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, all, all) = pressureData(vLevel)(stagCore(vLevel).ubound(0) - 1, all, all);
        }
#endif
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (not inputParams.yPer) {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION AT FRONT AND BACK WALLS
        if (zeroBC) {
            if (mesh.rankData.yRank == 0) {
                pressureData(vLevel)(all, -1, all) = -pressureData(vLevel)(all, 1, all);
            }

            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                pressureData(vLevel)(all, stagCore(vLevel).ubound(1) + 1, all) = -pressureData(vLevel)(all, stagCore(vLevel).ubound(1) - 1, all);
            }
        } else {
            if (mesh.rankData.yRank == 0) {
                pressureData(vLevel)(all, -1, all) = 2.0*yWall(all, all) - pressureData(vLevel)(all, 1, all);
            }

            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                pressureData(vLevel)(all, stagCore(vLevel).ubound(1) + 1, all) = 2.0*yWall(all, all) - pressureData(vLevel)(all, stagCore(vLevel).ubound(1) - 1, all);
            }
        }
#else
        // NEUMANN BOUNDARY CONDITION AT FRONT AND BACK WALLS
        if (mesh.rankData.yRank == 0) {
            pressureData(vLevel)(all, -1, all) = pressureData(vLevel)(all, 1, all);
        }

        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
            pressureData(vLevel)(all, stagCore(vLevel).ubound(1) + 1, all) = pressureData(vLevel)(all, stagCore(vLevel).ubound(1) - 1, all);
        }
#endif
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (inputParams.zPer) {
        // PERIODIC BOUNDARY CONDITION AT BOTTOM WALL
        pressureData(vLevel)(all, all, -1) = pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) - 1);

        // PERIODIC BOUNDARY CONDITION AT TOP WALL
        pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) + 1) = pressureData(vLevel)(all, all, 1);

    } else {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION AT BOTTOM AND TOP WALLS
        if (zeroBC) {
            pressureData(vLevel)(all, all, -1) = -pressureData(vLevel)(all, all, 1);

            pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) + 1) = -pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) - 1);
        } else {
            pressureData(vLevel)(all, all, -1) = 2.0*zWall(all, all) - pressureData(vLevel)(all, all, 1);

            pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) + 1) = 2.0*zWall(all, all) - pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) - 1);
        }
#else
        // NEUMANN BOUNDARY CONDITION AT BOTTOM AND TOP WALLS
        pressureData(vLevel)(all, all, -1) = pressureData(vLevel)(all, all, 1);

        pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) + 1) = pressureData(vLevel)(all, all, stagCore(vLevel).ubound(2) - 1);
#endif
    }
}


void multigrid_d3::updatePads(blitz::Array<blitz::Array<real, 3>, 1> &data) {
    recvRequest = MPI_REQUEST_NULL;

    // TRANSFER DATA FROM NEIGHBOURING CELL TO IMPOSE SUB-DOMAIN BOUNDARY CONDITIONS
    MPI_Irecv(&(data(vLevel)(mgRecvLft(vLevel))), 1, xMGArray(vLevel), mesh.rankData.faceRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(&(data(vLevel)(mgRecvRgt(vLevel))), 1, xMGArray(vLevel), mesh.rankData.faceRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));
    MPI_Irecv(&(data(vLevel)(mgRecvFrn(vLevel))), 1, yMGArray(vLevel), mesh.rankData.faceRanks(2), 3, MPI_COMM_WORLD, &recvRequest(2));
    MPI_Irecv(&(data(vLevel)(mgRecvBak(vLevel))), 1, yMGArray(vLevel), mesh.rankData.faceRanks(3), 4, MPI_COMM_WORLD, &recvRequest(3));

    MPI_Send(&(data(vLevel)(mgSendLft(vLevel))), 1, xMGArray(vLevel), mesh.rankData.faceRanks(0), 2, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendRgt(vLevel))), 1, xMGArray(vLevel), mesh.rankData.faceRanks(1), 1, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendFrn(vLevel))), 1, yMGArray(vLevel), mesh.rankData.faceRanks(2), 4, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendBak(vLevel))), 1, yMGArray(vLevel), mesh.rankData.faceRanks(3), 3, MPI_COMM_WORLD);

    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());

    // EVER SINCE THE DEVELOPMENT OF THE SOLVER STARTED IN 2018, THE ABOVE 4 SEND/IRECV CALLS WERE
    // SUFFICIENT FOR UPDATING THE PADS IN THE POISSON SOLVER. HOWEVER, IN JAN 2021, IT WAS FINALLY
    // DISCOVERED THAT THESE 4 CALLS WERE *NOT* ENOUGH, AND THAT THE VALUE AT THE SHARED POINTS
    // BETWEEN MPI SUB-DOMAIN BOUNDARIES NEEDED TO BE AVERAGED FROM BOTH THE SUB-DOMAINS.
    // THIS BUG MADE ITS APPEARANCE ONLY WHEN SOLVING CONVECTION PROBLEMS ON NON-UNIFORM GRIDS.
    MPI_Irecv(&(data(vLevel)(mgRecvLftFrn(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(&(data(vLevel)(mgRecvLftBak(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));
    MPI_Irecv(&(data(vLevel)(mgRecvRgtFrn(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(2), 3, MPI_COMM_WORLD, &recvRequest(2));
    MPI_Irecv(&(data(vLevel)(mgRecvRgtBak(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(3), 4, MPI_COMM_WORLD, &recvRequest(3));

    MPI_Send(&(data(vLevel)(mgSendLftFrn(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(0), 4, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendLftBak(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(1), 3, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendRgtFrn(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(2), 2, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendRgtBak(vLevel))), 1, zMGArray(vLevel), mesh.rankData.edgeRanks(3), 1, MPI_COMM_WORLD);

    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());
}


real multigrid_d3::testProlong() {
    vLevel = 0;

    // Fill the residualData array with correct values expected after prolongation
    residualData(vLevel) = 0.0;
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int j = 0; j <= yEnd(vLevel); ++j) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                residualData(vLevel)(i, j, k) = (mesh.rankData.rank + 1)*1000 + i*100 + j*10 + k;
            }
        }
    }

    // After going one level down the V-Cycle, populate the pressureData array with values at the corresponding stride
    vLevel += 1;
    pressureData(vLevel) = 0.0;
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int j = 0; j <= yEnd(vLevel); ++j) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                pressureData(vLevel)(i, j, k) = (mesh.rankData.rank + 1)*1000 + i*100 + j*10 + k;
            }
        }
    }

    // Perform prolongation
    prolong();

    pressureData(vLevel) -= residualData(vLevel - 1);

    return blitz::max(fabs(pressureData(vLevel)));
}


real multigrid_d3::testTransfer() {
    real maxVal = 0.0;

    vLevel = 0;

    pressureData(vLevel) = 0.0;
    residualData(vLevel) = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int j = 0; j <= yEnd(vLevel); ++j) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                pressureData(vLevel)(i, j, k) = (mesh.rankData.rank + 1)*1000 + i*100 + j*10 + k;
                residualData(vLevel)(i, j, k) = pressureData(vLevel)(i, j, k);
            }
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int j = 0; j <= yEnd(n); ++j) {
            for (int k = 0; k <= zEnd(n); ++k) {
                residualData(n)(-1, j, k) = (mesh.rankData.faceRanks(0) + 1)*1000 + (xEnd(n) - 1)*100 + j*10 + k;
                residualData(n)(xEnd(n) + 1, j, k) = (mesh.rankData.faceRanks(1) + 1)*1000 + 100 + j*10 + k;
            }
        }
    }

    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int i = 0; i <= xEnd(n); ++i) {
            for (int k = 0; k <= zEnd(n); ++k) {
                residualData(n)(i, -1, k) = (mesh.rankData.faceRanks(2) + 1)*1000 + i*100 + (yEnd(n) - 1)*10 + k;
                residualData(n)(i, yEnd(n) + 1, k) = (mesh.rankData.faceRanks(3) + 1)*1000 + i*100 + 10 + k;
            }
        }
    }

    for (int n=0; n<=inputParams.vcDepth; n++) {
        updatePads(pressureData);
        vLevel += 1;
    }

    pressureData(vLevel) -= residualData(vLevel);

    for (int i = pressureData(vLevel).lbound(0); i <= pressureData(vLevel).ubound(0); i += 1) {
        for (int j = pressureData(vLevel).lbound(1); j <= pressureData(vLevel).ubound(1); j += 1) {
            for (int k = 0; k <= zEnd(vLevel); k += 1) {
                if (abs(pressureData(vLevel)(i, j, k)) > maxVal) {
                    maxVal = abs(pressureData(vLevel)(i, j, k));
                }
            }
        }
    }

    return maxVal;
}


real multigrid_d3::testPeriodic() {
    real xCoord = 0.0;
    real yCoord = 0.0;
    real zCoord = 0.0;

    pressureData(0) = 0.0;
    residualData(0) = 0.0;

    for (int i = 0; i <= xEnd(0); ++i) {
        for (int j = 0; j <= yEnd(0); ++j) {
            for (int k = 0; k <= zEnd(0); ++k) {
                pressureData(0)(i, j, k) = sin(2.0*M_PI*mesh.x(i)/mesh.xLen)*
                                           cos(2.0*M_PI*mesh.y(j)/mesh.yLen)*
                                           cos(2.0*M_PI*mesh.z(k)/mesh.zLen);
                residualData(0)(i, j, k) = pressureData(0)(i, j, k);
            }
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int j = 0; j <= yEnd(n); ++j) {
            for (int k = 0; k <= zEnd(n); ++k) {
                xCoord = mesh.x(0) - (mesh.x(strideValues(n)) - mesh.x(0));
                residualData(n)(-1, j, k) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                            cos(2.0*M_PI*mesh.y(j)/mesh.yLen)*
                                            cos(2.0*M_PI*mesh.z(k)/mesh.zLen);

                xCoord = mesh.x(xEnd(0)) + (mesh.x(xEnd(0)) - mesh.x(xEnd(0) - strideValues(n)));
                residualData(n)(xEnd(n) + 1, j, k) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                     cos(2.0*M_PI*mesh.y(j)/mesh.yLen)*
                                                     cos(2.0*M_PI*mesh.z(k)/mesh.zLen);
            }
        }
    }

    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int i = 0; i <= xEnd(n); ++i) {
            for (int k = 0; k <= zEnd(n); ++k) {
                yCoord = mesh.y(0) - (mesh.y(strideValues(n)) - mesh.y(0));
                residualData(n)(i, -1, k) = sin(2.0*M_PI*mesh.x(i)/mesh.xLen)*
                                            cos(2.0*M_PI*yCoord/mesh.yLen)*
                                            cos(2.0*M_PI*mesh.z(k)/mesh.zLen);

                yCoord = mesh.y(yEnd(0)) + (mesh.y(yEnd(0)) - mesh.y(yEnd(0) - strideValues(n)));
                residualData(n)(i, yEnd(n) + 1, k) = sin(2.0*M_PI*mesh.x(i)/mesh.xLen)*
                                                     cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                     cos(2.0*M_PI*mesh.z(k)/mesh.zLen);
            }
        }
    }

    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int i = 0; i <= xEnd(n); ++i) {
            for (int j = 0; j <= yEnd(n); ++j) {
                zCoord = mesh.z(0) - (mesh.z(strideValues(n)) - mesh.z(0));
                residualData(n)(i, j, -1) = sin(2.0*M_PI*mesh.x(i)/mesh.xLen)*
                                            cos(2.0*M_PI*mesh.y(j)/mesh.yLen)*
                                            cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.z(zEnd(0)) + (mesh.z(zEnd(0)) - mesh.z(zEnd(0) - strideValues(n)));
                residualData(n)(i, j, zEnd(n) + 1) = sin(2.0*M_PI*mesh.x(i)/mesh.xLen)*
                                                     cos(2.0*M_PI*mesh.y(j)/mesh.yLen)*
                                                     cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    vLevel = 0;
    for (int n=0; n<=inputParams.vcDepth; n++) {
        imposeBC();
        vLevel += 1;
    }

    pressureData(vLevel) -= residualData(vLevel);

    return blitz::max(fabs(pressureData(vLevel)));
}
