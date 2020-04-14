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
/*! \file poisson.cc
 *
 *  \brief Definitions for functions of class poisson
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
 * \brief   Constructor of the base poisson class
 *
 *          The short base constructor of the poisson class merely assigns the const references to the grid and parser
 *          class instances being used in the solver.
 *          Moreover, it resizes and populates a local array of multi-grid sizes as used in the grid class.
 *          An array of strides to be used at different V-cycle levels is also generated and stored.
 *          Finally, the maximum allowable number of iterations for the Jacobi iterative solver being used at the
 *          coarsest mesh is set as \f$ N_{max} = N_x \times N_y \times N_z \f$, where \f$N_x\f$, \f$N_y\f$ and \f$N_z\f$
 *          are the number of grid points in the collocated grid at the local sub-domains along x, y and z directions
 *          respectively.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
poisson::poisson(const grid &mesh, const parser &solParam): mesh(mesh), inputParams(solParam) {
    int maxIndex = 15;

    mgSizeArray.resize(maxIndex);
    for (int i=0; i < maxIndex; i++) {
        mgSizeArray(i) = int(pow(2, i)) + 1;
    }

    mgSizeArray(0) = 1;

    strideValues.resize(inputParams.vcDepth + 1);
    for (int i=0; i<=inputParams.vcDepth; i++) {
        strideValues(i) = int(pow(2, i));
    }

    vLevel = 0;
    maxCount = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);

#ifdef TIME_RUN
    smothTimeComp = 0.0;
    smothTimeTran = 0.0;
#endif
}


/**
 ********************************************************************************************************************************************
 * \brief   The core, publicly accessible function of poisson to compute the solution for the Poisson equation
 *
 *          The function calls the V-cycle as many times as set by the user.
 *          Before doing so, the input data is transferred into the data-structures used by the poisson class to
 *          perform restrictions and prolongations without copying.
 *          Finally, the computed solution is transferred back from the internal data-structures back into the
 *          scalar field supplied by the calling function.
 *
 * \param   inFn is a pointer to the plain scalar field (cell-centered) into which the computed soltuion must be transferred
 * \param   rhs is a const reference to the plain scalar field (cell-centered) which contains the RHS for the Poisson equation to solve
 ********************************************************************************************************************************************
 */
void poisson::mgSolve(plainsf &inFn, const plainsf &rhs) {
#ifdef TEST_POISSON
    // The below value can be chosen to compute the residual in one of the two ways:
    // 0 = Mean Absolute Error (MAE)
    // 1 = Root Mean Squared Error (RMSE)
    int normType = 1;
    real mgResidual;
#endif

    pressureData = 0.0;
    residualData = 0.0;
    inputRHSData = 0.0;

    // TRANSFER DATA FROM THE INPUT SCALAR FIELDS INTO THE DATA-STRUCTURES USED BY poisson
    inputRHSData(stagCore) = rhs.F(stagCore);
    pressureData(stagCore) = inFn.F(stagCore);

#ifndef TEST_POISSON
    // TO MAKE THE PROBLEM WELL-POSED (WHEN USING NEUMANN BC ONLY), SUBTRACT THE MEAN OF THE RHS FROM THE RHS
    real localMean = blitz::mean(inputRHSData);
    real globalAvg = 0.0;

    MPI_Allreduce(&localMean, &globalAvg, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
    globalAvg /= mesh.rankData.nProc;

    inputRHSData -= globalAvg;
#endif

    // PERFORM V-CYCLES AS MANY TIMES AS REQUIRED
    for (int i=0; i<inputParams.vcCount; i++) {
        smoothedPres = 0.0;

        vCycle();

#ifdef TEST_POISSON
        mgResidual = computeError(normType);
        if (mesh.rankData.rank == 0) std::cout << "Residual after V Cycle is " << mgResidual << std::endl;
#endif
    }

    // RETURN CALCULATED PRESSURE DATA
    inFn.F = pressureData(blitz::RectDomain<3>(inFn.F.lbound(), inFn.F.ubound()));
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to perform one loop of V-cycle
 *
 *          The V-cycle of restrictions, prolongations and smoothings are performed within this function.
 *          First the input data contained in \ref pressureData is smoothed, after which the residual is computed and stored
 *          in the \ref residualData array.
 *          The restrictions, smoothing, and prolongations are performed on these two arrays subsequently.
 ********************************************************************************************************************************************
 */
void poisson::vCycle() {
    /*
     * OUTLINE OF THE MULTI-GRID V-CYCLE
     * 1)  Start at finest grid, perform N=2 Gauss-Siedel pre-smoothing iterations to solve for the solution Ax=b,
     * 2)  Compute the residual r=b-Ax and restrict it to a coarser level,
     * 3)  Perform N=2 Gauss-Siedel pre-smoothing iterations to solve for the error: Ae=r
     * 4)  Repeat steps 2-3 until you reach the coarsest grid level,
     * 5)  Perform N=2+2 (pre+post) Gauss-Siedel smoothing iterations to solve for the error 'e',
     * 6)  Prolong the error 'e' to the next finer level.
     * 7)  Perform N=2 post-smoothing iterations, 
     * 8)  Repeat steps 6-7 until the finest grid is reached,
     * 9)  Add error 'e' to the solution 'x' and perform N=2 post-smoothing iterations.
     * 10) End of one V-cycle, and check for convergence by computing the normalized residual: r_normalized = ||b-Ax||/||b||. 
     */

    vLevel = 0;

    //std::cout << "Pre-smoothing" << std::endl;
    // Step 1) Pre-smoothing iterations of Ax = b
    swap(inputRHSData, residualData);
    smooth(inputParams.preSmooth);
    swap(residualData, inputRHSData);
    // After above 3 lines, pressureData has the pre-smoothed values of pressure, inputRHSData has original RHS data, and residualData = 0.0

    // Step 2) Compute the residual r = b - Ax
    computeResidual();

    //std::cout << residualData(blitz::Range::all(), 5, 5) << std::endl;
    //MPI_Finalize();
    //exit(0);
    // Shift pressureData into smoothedPres
    swap(smoothedPres, pressureData);
    // Now pressureData = 0.0 (since smoothedPres was 0.0), and smoothedPres has the pre-smoothed values of pressure
    // So smoothedPres contains smoothed x, residualData holds r, and pressureData is ready to hold e.

    // RESTRICTION OPERATIONS
    for (int i=0; i<inputParams.vcDepth; i++) {
        coarsen();
    //std::cout << "Coarsening-smoothing" << i << std::endl;
        // Step 3) Perform pre-smoothing iterations to solve for the error: Ae = r
        smooth(inputParams.restrictSmooth[i]);
    }
    // Step 4) Repeat steps 2-3 until you reach the coarsest grid level,

    //std::cout << "Final-smoothing" << std::endl;
    // Step 5) Perform pre+post smoothing iterations to solve for the error 'e',
    smooth(inputParams.restrictSmooth[inputParams.vcDepth] + inputParams.prolongSmooth[inputParams.vcDepth]);

    // PROLONGATION OPERATIONS BACK TO FINE MESH
    for (int i=0; i<inputParams.vcDepth; i++) {
        // Step 6) Prolong the error 'e' to the next finer level.
        prolong();
    //std::cout << "Prolong-smoothing" << i << std::endl;
        smooth(inputParams.prolongSmooth[i]);
    }

    // Step 9) Add error 'e' to the solution 'x' and perform post-smoothing iterations.
    pressureData += smoothedPres;

    // POST-SMOOTHING
    swap(inputRHSData, residualData);
    //std::cout << "Post-smoothing" << std::endl;
    smooth(inputParams.postSmooth);
    swap(residualData, inputRHSData);

    //MPI_Finalize();
    //exit(0);
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the arrays used in multi-grid
 *
 *          The memory required for various arrays in multi-grid solver are pre-allocated through this function.
 *          The function is called from within the constructor to perform this allocation once and for all.
 *          The arrays are initialized to 0.
 ********************************************************************************************************************************************
 */
void poisson::initializeArrays() {
    pressureData.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    pressureData.reindexSelf(stagFull.lbound());
    pressureData = 0.0;

    smoothedPres.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    smoothedPres.reindexSelf(stagFull.lbound());
    smoothedPres = 0.0;

    iteratorTemp.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    iteratorTemp.reindexSelf(stagFull.lbound());
    iteratorTemp = 0.0;

    inputRHSData.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    inputRHSData.reindexSelf(stagFull.lbound());
    inputRHSData = 0.0;

    residualData.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    residualData.reindexSelf(stagFull.lbound());
    residualData = 0.0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the Range objects for accessing mesh derivatives in transformed plane
 *
 *          The Range objects defined here are used for reading the values of grid metrics at all the V-cycle levels.
 *          Since these values are required at all the grid levels, there are \ref parser#vcDepth "vcDepth" + 1 number of Range objects.
 ********************************************************************************************************************************************
 */
void poisson::initMeshRanges() {
    xMeshRange.resize(inputParams.vcDepth + 1);
#ifndef PLANAR
    yMeshRange.resize(inputParams.vcDepth + 1);
#endif
    zMeshRange.resize(inputParams.vcDepth + 1);

    // Range OBJECTS WITH STRIDE TO ACCESS DIFFERENT POINTS OF THE SAME ARRAY AT DIFFERENT MULTI-GRID LEVELS
    for(int i=0; i<=inputParams.vcDepth; i++) {
        xMeshRange(i) = blitz::Range(stagCore.lbound(0), stagCore.ubound(0), strideValues(i));
#ifndef PLANAR
        yMeshRange(i) = blitz::Range(stagCore.lbound(1), stagCore.ubound(1), strideValues(i));
#endif
        zMeshRange(i) = blitz::Range(stagCore.lbound(2), stagCore.ubound(2), strideValues(i));
    }

    // SET THE LIMTS FOR ARRAY LOOPS IN smooth FUNCTION, AND A FEW OTHER PLACES
    // WARNING: THESE VARIABLES HAVE SO FAR BEEN IMPLEMENTED ONLY IN smooth AND vCycle.
    // THE TEST FUNCTIONS HAVE NOT YET BEEN UPDATED WITH THESE
    xStr = stagCore.lbound(0);
    xEnd = stagCore.ubound(0);

#ifndef PLANAR
    yStr = stagCore.lbound(1);
    yEnd = stagCore.ubound(1);
#endif

    zStr = stagCore.lbound(2);
    zEnd = stagCore.ubound(2);
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the RectDomain variables for all future references throughout the poisson solver
 *
 *          The function sets the core and full domain staggered grid sizes for all the sub-domains.
 ********************************************************************************************************************************************
 */
void poisson::setStagBounds() {
    blitz::TinyVector<int, 3> loBound, upBound;

    // LOWER BOUND AND UPPER BOUND OF STAGGERED CORE - USED TO CONSTRUCT THE CORE SLICE
    loBound = 0, 0, 0;
#ifdef PLANAR
    upBound = mgSizeArray(localSizeIndex(0)) - 1, 0, mgSizeArray(localSizeIndex(2)) - 1;
#else
    upBound = mgSizeArray(localSizeIndex(0)) - 1, mgSizeArray(localSizeIndex(1)) - 1, mgSizeArray(localSizeIndex(2)) - 1;
#endif
    stagCore = blitz::RectDomain<3>(loBound, upBound);

    // LOWER BOUND AND UPPER BOUND OF STAGGERED FULL SUB-DOMAIN - USED TO CONSTRUCT THE FULL SUB-DOMAIN SLICE
#ifdef PLANAR
    loBound = -strideValues(inputParams.vcDepth), -1, -strideValues(inputParams.vcDepth);
#else
    loBound = -strideValues(inputParams.vcDepth), -strideValues(inputParams.vcDepth), -strideValues(inputParams.vcDepth);
#endif
    upBound = stagCore.ubound() - loBound;
    stagFull = blitz::RectDomain<3>(loBound, upBound);
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the local size indices of sub-domains after MPI domain decomposition
 *
 *          For the multi-grid solver, the number of grid nodes must be \f$ 2^N + 1 \f$ to perform V-Cycles
 *          The domain decomposition is also done in such a way that each sub-domain will also have \f$ 2^M + 1 \f$ points.
 *          As a result, the number of processors in each direction must be a power of 2.
 *          In this case, the the value \f$ M \f$ of \f$ 2^M + 1 \f$ can be computed as \f$ M = N - log_2(np) \f$
 *          where \f$ np \f$ is the number of processors along the direction under consideration.
 ********************************************************************************************************************************************
 */
void poisson::setLocalSizeIndex() {
#ifdef PLANAR
    localSizeIndex = blitz::TinyVector<int, 3>(mesh.sizeIndex(0) - int(log2(inputParams.npX)),
                                               mesh.sizeIndex(1),
                                               mesh.sizeIndex(2));
#else
    localSizeIndex = blitz::TinyVector<int, 3>(mesh.sizeIndex(0) - int(log2(inputParams.npX)),
                                               mesh.sizeIndex(1) - int(log2(inputParams.npY)),
                                               mesh.sizeIndex(2));
#endif
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the coefficients used for calculating laplacian and in smoothing
 *
 *          The function assigns values to the variables \ref hx, \ref hy etc.
 *          These coefficients are repeatedly used at many places in the Poisson solver.
 ********************************************************************************************************************************************
 */
void poisson::setCoefficients() {
    hx.resize(inputParams.vcDepth + 1);
#ifndef PLANAR
    hy.resize(inputParams.vcDepth + 1);
#endif
    hz.resize(inputParams.vcDepth + 1);

#ifdef PLANAR
    hx2.resize(inputParams.vcDepth + 1);
    hz2.resize(inputParams.vcDepth + 1);
#else
    hxhy.resize(inputParams.vcDepth + 1);
    hyhz.resize(inputParams.vcDepth + 1);
#endif
    hzhx.resize(inputParams.vcDepth + 1);

#ifndef PLANAR
    hxhyhz.resize(inputParams.vcDepth + 1);
#endif

    for(int i=0; i<=inputParams.vcDepth; i++) {
        hx(i) = strideValues(i)*mesh.dXi;
#ifndef PLANAR
        hy(i) = strideValues(i)*mesh.dEt;
#endif
        hz(i) = strideValues(i)*mesh.dZt;

#ifdef PLANAR
        hx2(i) = pow(strideValues(i)*mesh.dXi, 2.0);
        hz2(i) = pow(strideValues(i)*mesh.dZt, 2.0);
#else
        hxhy(i) = pow(strideValues(i), 4.0)*pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0);
        hyhz(i) = pow(strideValues(i), 4.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
#endif
        hzhx(i) = pow(strideValues(i), 4.0)*pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);

#ifndef PLANAR
        hxhyhz(i) = pow(strideValues(i), 6.0)*pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
#endif
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to copy the staggered grid derivatives from the grid class to local arrays
 *
 *          Though the grid derivatives in the grid class can be read and accessed, they cannot be used directly
 *          along with the arrays defined in the poisson class as the local arrays have wide pads on both sides.
 *          Therefore, correspondingly wide arrays for grid derivatives are used, into which the staggered grid
 *          derivatives from the grid class are written and stored.
 *          This function serves this purpose of copying the grid derivatives.
 ********************************************************************************************************************************************
 */
void poisson::copyStaggrDerivs() {
    xixx.resize(stagFull.ubound(0) - stagFull.lbound(0) + 1);
    xixx.reindexSelf(stagFull.lbound(0));
    xixx = 0.0;
    xixx(blitz::Range(0, stagCore.ubound(0), 1)) = mesh.xixxStaggr(blitz::Range(0, stagCore.ubound(0), 1));

    xix2.resize(stagFull.ubound(0) - stagFull.lbound(0) + 1);
    xix2.reindexSelf(stagFull.lbound(0));
    xix2 = 0.0;
    xix2(blitz::Range(0, stagCore.ubound(0), 1)) = mesh.xix2Staggr(blitz::Range(0, stagCore.ubound(0), 1));

#ifndef PLANAR
    etyy.resize(stagFull.ubound(1) - stagFull.lbound(1) + 1);
    etyy.reindexSelf(stagFull.lbound(1));
    etyy = 0.0;
    etyy(blitz::Range(0, stagCore.ubound(1), 1)) = mesh.etyyStaggr(blitz::Range(0, stagCore.ubound(1), 1));

    ety2.resize(stagFull.ubound(1) - stagFull.lbound(1) + 1);
    ety2.reindexSelf(stagFull.lbound(1));
    ety2 = 0.0;
    ety2(blitz::Range(0, stagCore.ubound(1), 1)) = mesh.ety2Staggr(blitz::Range(0, stagCore.ubound(1), 1));
#endif

    ztzz.resize(stagFull.ubound(2) - stagFull.lbound(2) + 1);
    ztzz.reindexSelf(stagFull.lbound(2));
    ztzz = 0.0;
    ztzz(blitz::Range(0, stagCore.ubound(2), 1)) = mesh.ztzzStaggr(blitz::Range(0, stagCore.ubound(2), 1));

    ztz2.resize(stagFull.ubound(2) - stagFull.lbound(2) + 1);
    ztz2.reindexSelf(stagFull.lbound(2));
    ztz2 = 0.0;
    ztz2(blitz::Range(0, stagCore.ubound(2), 1)) = mesh.ztz2Staggr(blitz::Range(0, stagCore.ubound(2), 1));
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to coarsen the grid down the levels of the V-Cycle
 *
 *          Coarsening reduces the number of points in the grid by averaging values at two adjacent nodes onto an intermediate point between them
 *          As a result, the number of points in the domain decreases from \f$ 2^{N+1} + 1 \f$ at the input level to \f$ 2^N + 1 \f$.
 *          The vLevel variable is accordingly incremented by 1 to reflect this descent by one step down the V-Cycle.
 ********************************************************************************************************************************************
 */
void poisson::coarsen() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to perform prolongation on the array being solved
 *
 *          Prolongation makes the grid finer by averaging values at two adjacent nodes onto an intermediate point between them
 *          As a result, the number of points in the domain increases from \f$ 2^N + 1 \f$ at the input level to \f$ 2^{N+1} + 1 \f$.
 *          The vLevel variable is accordingly reduced by 1 to reflect this ascent by one step up the V-Cycle.
 ********************************************************************************************************************************************
 */
void poisson::prolong() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the residual at the start of each V-Cycle
 *
 *          The Poisson solver solves for the residual r = b - Ax
 *          This function computes this residual by calculating the Laplacian of the pressure field and
 *          subtracting it from the RHS of Poisson equation.
 *
 ********************************************************************************************************************************************
 */
void poisson::computeResidual() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to perform smoothing operation on the input array
 *
 *          The smoothing operation is always performed on the data contained in the array \ref pressureData.
 *          The array \ref iteratorTemp is used to store the temporary data and it is continuously swapped with the
 *          \ref pressureData array at every iteration.
 *          This operation can be performed at any level of the V-cycle.
 *
 * \param   smoothCount is the integer value of the number of smoothing iterations to be performed
 ********************************************************************************************************************************************
 */
void poisson::smooth(const int smoothCount) { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the error at the end of each V-Cycle
 *
 *          To check for convergence, the residual must be computed throughout the domain.
 *          This function offers multiple ways to compute the residual (global maximum, rms, mean, etc.)
 *
 * \param   normOrder is the integer value of the order of norm used to calculate the residual.
 ********************************************************************************************************************************************
 */
real poisson::computeError(const int normOrder) {  return 0.0; };


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions of Poisson solver at different levels of the V-cycle
 *
 *          This function is called mainly during smoothing operations to impose the boundary conditions for the
 *          Poisson equation.
 *          The sub-domains close to the wall will have the Neumann boundary condition on pressure imposeed at the walls.
 *          Meanwhile at the interior boundaries at the inter-processor sub-domains, data is transferred from the neighbouring cells
 *          by calling the \ref updatePads function.
 ********************************************************************************************************************************************
 */
void poisson::imposeBC() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to update the pad points of the local sub-domains at different levels of the V-cycle
 *
 *          This function is called mainly during smoothing operations by the \ref imposeBC function.
 *          At the interior boundaries at the inter-processor sub-domains, data is transferred from the neighbouring cells
 *          using a combination of MPI_Irecv and MPI_Send functions.
 ********************************************************************************************************************************************
 */
void poisson::updatePads() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to create the MPI sub-array data types necessary to transfer data across sub-domains
 *
 *          The inter-domain boundaries of all the sub-domains at different V-cycle levels need data to be transfered at
 *          with different mesh strides.
 *          The number of sub-arrays along each edge/face of the sub-domains are equal to the number of V-cycle levels.
 *          Since this data transfer has to take place at all the mesh levels including the finest mesh, there will be
 *          vcDepth + 1 elements.
 ********************************************************************************************************************************************
 */
void poisson::createMGSubArrays() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether strided data transfer is performing as expected
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls updatePads at different vLevels and checks is the data is being transferred along x and y directions
 *          This done by printing the contents of the arrays for visual inspection for now.
 ********************************************************************************************************************************************
 */
real poisson::testTransfer() { return 0; };

/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether prolongation operations are interpolating correctly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls prolong function at a lower vLevel and checks if the data is being interpolated correctly at higher vLevel
 *          This done by returning the average deviation from correct values as a real value
 ********************************************************************************************************************************************
 */
real poisson::testProlong() { return 0; };

/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether periodic BC is being implemented properly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls imposeBC function at different vLevels and checks if the correct values of the functions are imposed at boundaries
 *          This done by printing the contents of the arrays for visual inspection for now.
 ********************************************************************************************************************************************
 */
real poisson::testPeriodic() { return 0; };

poisson::~poisson() {
#ifdef TIME_RUN
    if (mesh.rankData.rank == 0) {
        std::cout << std::left << std::setw(50) << "Time taken in computation within smooth: "           << std::fixed << std::setprecision(6) << smothTimeComp << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in data-transfer within smooth: "         << std::fixed << std::setprecision(6) << smothTimeTran << std::endl;
    }
#endif
};
