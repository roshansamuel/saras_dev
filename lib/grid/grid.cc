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
/*! \file grid.cc
 *
 *  \brief Definitions for functions of class grid
 *  \sa grid.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "grid.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the grid class
 *
 *          The class constructor initializes the mesh for computational problem.
 *          The pad widths, global grid limits in the full domain, local grid limits in the MPI decomposed sub-domains,
 *          grid spacings, domain lengths, etc., along each direction are set.
 *          Appropriate stretching functions are chosen according to user preferences and their corresponding grid
 *          transformation derivatives are also computed and stored.
 *
 * \param   solParam is a const reference to the global data contained in the parser class
 * \param   parallelData is a reference to the global data contained in the parallel class
 ********************************************************************************************************************************************
 */
grid::grid(const parser &solParam, parallel &parallelData): inputParams(solParam),
                                                            rankData(parallelData) {
    /** Depending on the finite-difference scheme chosen for calculating derivatives, set the \ref padWidths along all directions. */
    if (inputParams.dScheme == 1) {
        padWidths = 1, 1, 1;
    } else if (inputParams.dScheme == 2) {
        padWidths = 2, 2, 2;
    } else {
        if (rankData.rank == 0) {
            std::cout << "Undefined finite differencing scheme in YAML file. ABORTING" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    // THE ARRAY sizeArray HAS ELEMENTS [1, 3, 5, 9, 17, 33 ..... ] - STAGGERED GRID SIZE
    makeSizeArray();

    sizeIndex = inputParams.xInd, inputParams.yInd, inputParams.zInd;
    globalSize = sizeArray(sizeIndex(0)), sizeArray(sizeIndex(1)), sizeArray(sizeIndex(2));
#ifdef PLANAR
    totalPoints = globalSize(0)*globalSize(2);
#else
    totalPoints = globalSize(0)*globalSize(1)*globalSize(2);
#endif

    xLen = inputParams.Lx;
    yLen = inputParams.Ly;
    zLen = inputParams.Lz;

    thBeta = inputParams.betaX, inputParams.betaY, inputParams.betaZ;

    dXi = 1.0/real(globalSize(0) - 1);
    dEt = 1.0/real(globalSize(1) - 1);
    dZt = 1.0/real(globalSize(2) - 1);

#ifdef PLANAR
    // IS IT OKAY TO SET BELOW VALUE AS 1 EVEN WHEN inputParams.dScheme IS NOT 1?
    padWidths(1) = 1;
    yLen = 1.0;
    dEt = 1.0;
#endif

    // COMPUTE THE LOCAL ARRAY SIZES, coreSize, START AND END INDICES, subarrayStarts AND subarrayEnds
    computeGlobalLimits();

    // SET THE TinyVector AND RectDomain VARIABLES BASED ON VALUES COMPUTED IN computeGlobalLimits, FOR RESIZING ALL LOCAL GRIDS
    setDomainSizes();

    // RESIZE GRID USING THE VARIABLES CONSTRUCTED ABOVE IN setDomainSizes
    resizeGrid();

    // GENERATE THE GLOBAL TRANSFORMED GRID
    globalXiEtaZeta();

    // SET LOCAL TRANSFORMED GRID AS SLICES FROM THE GLOBAL TRANSFORMED GRID GENERATED ABOVE IN globalXiEtaZeta
    xi = xiGlo(blitz::Range(subarrayStarts(0) - padWidths(0), subarrayEnds(0) + padWidths(0)));
    et = etGlo(blitz::Range(subarrayStarts(1) - padWidths(1), subarrayEnds(1) + padWidths(1)));
    zt = ztGlo(blitz::Range(subarrayStarts(2) - padWidths(2), subarrayEnds(2) + padWidths(2)));

    // CREATE UNIFORM GRID WHICH IS DEFAULT ALONG ALL THREE DIRECTIONS
    createUniformGrid();

    // FLAG TO CHECK FOR GRID ANISOTROPY - FALSE BY DEFAULT UNLESS NON-UNIFORM GRID IS CREATED
    bool gridCheck = false;

    // DEPENDING ON THE USER-SET PARAMETERS, SWITCH TO TAN-HYP ALONG SELECTED DIRECTIONS
    if (inputParams.xGrid == 2) {
        createTanHypGrid(0);
        gridCheck = true;
    }

    if (inputParams.yGrid == 2) {
        createTanHypGrid(1);
        gridCheck = true;
    }

    if (inputParams.zGrid == 2) {
        createTanHypGrid(2);
        gridCheck = true;
    }

    syncGrid();

    if (gridCheck) checkAnisotropy();

    gatherGlobal();
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to resize and initialize the size array from which the dimensions of the grid will be determined
 *
 *          The size array will generate grid sizes according to \f$ 2^N + 2 \f$ to enable multigrid operations on the grid.
 *          The \ref parser#xInd "xInd", \ref parser#yInd "yInd" and \ref parser#zInd "zInd" parameters set by the users in
 *          parameters.yaml and read by the \ref parser class will be used to locate the grid size within the \ref sizeArray and
 *          generate grid accordingly.
 *
 *          Note that the grid sizes stored in \ref sizeArray correspond to the collocated grid.
 *          For multi-grid operations, the number of grid points necessary is \f$ 2^N + 1 \f$, which is 1 less than the values
 *          generated by this function.
 *          However, multi-grid is applied here to compute pressure correction and pressure is calculated on the staggered grid.
 *          Since there are \f$ N - 1 \f$ staggered grid points for \f$ N \f$ collocated points, the sizes become consistent.
 ********************************************************************************************************************************************
 */
void grid::makeSizeArray() {
    int maxIndex = 15;

    sizeArray.resize(maxIndex);
    for (int i=0; i < maxIndex; i++) {
        sizeArray(i) = int(pow(2, i)) + 1;
    }

    sizeArray(0) = 1;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the extent of local sub-domains in terms of the global index of the full domain
 *
 *          Depending on the number of processor divisions along each direction, the limits of the grid for each local
 *          sub-domain is set based on its \ref parallel#xRank "xRank" and \ref parallel#yRank "yRank".
 *          These limits are used to locate the local sub-domains within the full domain later.
 ********************************************************************************************************************************************
 */
void grid::computeGlobalLimits() {
    int xiSt, etSt, ztSt;
    int xiEn, etEn, ztEn;
    int localNx, localNy, localNz;

    // NUMBER OF STAGGERED POINTS IN EACH SUB-DOMAIN EXCLUDING PAD POINTS
    localNx = (globalSize(0) - 1)/rankData.npX + 1;
#ifndef PLANAR
    localNy = (globalSize(1) - 1)/rankData.npY + 1;
#else
    localNy = 1;
#endif
    localNz = (globalSize(2));

    // SETTING GLOBAL LIMITS
    // ADD ONE EXTRA POINT EACH AT FIRST AND LAST SUB-DOMAINS
    // FIRST SET THE LIMITS TO DEFAULT VALUES - THIS ELIMINATES AN EXTRA 'if' CONDITION
    // THEN SET LIMITS FOR LAST RANK IN EACH DIRECTION FIRST AND *FINALLY* SET LIMITS OF 0TH RANK
    // THIS IS NECESSARY TO AVOID ERRORS WHEN A PROCESSOR IS BOTH FIRST AND LAST RANK
    // THIS HAPPENS WHEN THERE ARE NO DIVISIONS ALONG AN AXIS AS ALONG Z-DIRECTION

    // ALONG XI-DIRECTION
    xiSt = rankData.xRank*(localNx - 1);
    xiEn = xiSt + localNx - 1;

    // ALONG ETA-DIRECTION
    etSt = rankData.yRank*(localNy - 1);
    etEn = etSt + localNy - 1;

    // ALONG ZETA-DIRECTION
    ztSt = 0;
    ztEn = ztSt + localNz - 1;

    coreSize = localNx, localNy, localNz;
    fullSize = coreSize + 2*padWidths;

    // SUB-ARRAY STARTS AND ENDS FOR *STAGGERED* GRID
    subarrayStarts = xiSt, etSt, ztSt;
    subarrayEnds = xiEn, etEn, ztEn;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set all the TinyVector and RectDomain variables for all future references throughout the solver
 *
 *          The function sets the core and full domain sizes for all the sub-domains after MPI decomposition.
 *          Additionally, the pad widths and starting indices of the sub-domains within the global domain are also set.
 ********************************************************************************************************************************************
 */
void grid::setDomainSizes() {
    blitz::TinyVector<int, 3> loBound, upBound;

    // LOWER BOUND AND UPPER BOUND OF CORE - USED TO CONSTRUCT THE CORE SLICE OF STAGGERED POINTS
    upBound = coreSize - 1;
    coreDomain = blitz::RectDomain<3>(loBound, upBound);

    // LOWER BOUND AND UPPER BOUND OF FULL SUB-DOMAIN - USED TO CONSTRUCT THE FULL SUB-DOMAIN SLICE
    upBound = coreSize + padWidths - 1;
    fullDomain = blitz::RectDomain<3>(loBound, upBound);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to resize and initialize the grid
 *
 *          The global collocated grid in transformed plane are resized according the global size of full domain.
 *          Then the local collocated grid in transformed plane is resized according to the limits defined in \ref computeGlobalLimits.
 *          Correspondingly, this function is called after the global limits have been set.
 *          After defining the transformed plane coordinates, both the staggered and collocated grids in physical plane are resized.
 *          Finally, the arrays for the grid derivative terms are also resized and initialized to 1.
 ********************************************************************************************************************************************
 */
void grid::resizeGrid() {
    // ALL ARRAYS MUST BE RESIZED AND REINDEXED TO LET THE NEGATIVE PADS HAVE NEGATIVE INDICES
    // THIS IS DONE IN A SINGLE STEP BY INITIALIZING THE ARRAYS WITH A blitz::Range OBJECT WHICH CONTAINS SIZE AND INDEXING INFORMATION
    blitz::Range xRange, yRange, zRange;

    // RANGE OF THE SUB-DOMAIN FOR STAGGERED AND COLLOCATED GRIDS: CONSTRUCTED FROM LOWER AND UPPER BOUNDS OF FULL SUB-DOMAIN
    xRange = blitz::Range(fullDomain.lbound(0), fullDomain.ubound(0));
    yRange = blitz::Range(fullDomain.lbound(1), fullDomain.ubound(1));
    zRange = blitz::Range(fullDomain.lbound(2), fullDomain.ubound(2));

    // LOCAL XI, ETA AND ZETA ARRAYS
    xi.resize(xRange);
    et.resize(yRange);
    zt.resize(zRange);

    // STAGGERED GRID POINTS AND THEIR METRICS
    x.resize(xRange);
    y.resize(yRange);
    z.resize(zRange);

    xi_x.resize(xRange);        xixx.resize(xRange);        xix2.resize(xRange);
    et_y.resize(yRange);        etyy.resize(yRange);        ety2.resize(yRange);
    zt_z.resize(zRange);        ztzz.resize(zRange);        ztz2.resize(zRange);

    x = 1.0;        y = 1.0;        z = 1.0;

    // BELOW ARE DEFAULT VALUES FOR A UNIFORM GRID OVER DOMAIN OF LENGTH 1.0
    // THESE VALUES ARE OVERWRITTEN AS PER GRID TYPE
    xi_x = 1.0;     xixx = 0.0;     xix2 = 1.0;
    et_y = 1.0;     etyy = 0.0;     ety2 = 1.0;
    zt_z = 1.0;     ztzz = 0.0;     ztz2 = 1.0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute global values of xi, eta and zeta in transformed plane
 *
 *          The function populates the \ref xiGlo, \ref etGlo and \ref ztGlo arrays from which the local values of
 *          \ref xi, \ref et and \ref zt in each sub-domain are obtained.
 *          These local values are obtained from the global grid according to the limits defined in \ref computeGlobalLimits.
 ********************************************************************************************************************************************
 */
void grid::globalXiEtaZeta() {
    xiGlo.resize(globalSize(0) + 2*padWidths(0));          xiGlo.reindexSelf(-padWidths(0));
    etGlo.resize(globalSize(1) + 2*padWidths(1));          etGlo.reindexSelf(-padWidths(1));
    ztGlo.resize(globalSize(2) + 2*padWidths(2));          ztGlo.reindexSelf(-padWidths(2));

    // ALONG XI-DIRECTION
    for (int i=-padWidths(0); i<globalSize(0)+padWidths(0); i++) {
        xiGlo(i) = real(i)*dXi;
    }

    // ALONG ETA-DIRECTION
    for (int i=-padWidths(1); i<globalSize(1)+padWidths(1); i++) {
        etGlo(i) = real(i)*dEt;
    }

    // ALONG ZETA-DIRECTION
    for (int i=-padWidths(2); i<globalSize(2)+padWidths(2); i++) {
        ztGlo(i) = real(i)*dZt;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to generate grid with uniform stretching
 *
 *          The local collocated grids, \ref xColloc, \ref yColloc, and \ref zColloc are equated to their corresponding
 *          transformed plane coordinates, \ref xi, \ref et, and \ref zt respectively.
 *          The corresponding grid derivative terms, \ref xi_xColloc, \ref xixxColloc, \ref ety2Colloc, etc are left
 *          unchanged from their initial value of 1.0, indicating that the grid is uniform.
 *
 *          Similarly, the staggered grids, \ref xStaggr, \ref yStaggr, and \ref zStaggr are also equated to the mid-point
 *          averaged values of the nodes in their corresponding transformed plane coordinates, \ref xi, \ref et, and \ref zt
 *          respectively.
 *          As before, the grid derivative terms for the staggered points are also left as 1.0.
 ********************************************************************************************************************************************
 */
void grid::createUniformGrid() {
    int i;

    // COLLOCATED X-GRID POINTS FROM UNIFORM XI-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(0); i < coreSize(0) + padWidths(0); i++) {
        x(i) = xLen*xi(i);
    }

#ifndef PLANAR
    // COLLOCATED Y-GRID POINTS FROM UNIFORM ETA-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(1); i < coreSize(1) + padWidths(1); i++) {
        y(i) = yLen*et(i);
    }
#endif

    // COLLOCATED Z-GRID POINTS FROM UNIFORM ZETA-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(2); i < coreSize(2) + padWidths(2); i++) {
        z(i) = zLen*zt(i);
    }

    xi_x = 1.0/xLen;
    xix2 = pow(xi_x, 2.0);

    et_y = 1.0/yLen;
    ety2 = pow(et_y, 2.0);

    zt_z = 1.0/zLen;
    ztz2 = pow(zt_z, 2.0);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to generate grid with tangent-hyperbolic stretching
 *
 *          The local collocated grids, \ref xColloc, \ref yColloc, and \ref zColloc are initialized from their corresponding
 *          transformed plane coordinates, \ref xi, \ref et, and \ref zt respectively using the tangent hyperbolic function.
 *          The corresponding grid derivative terms, \ref xi_xColloc, \ref xixxColloc, \ref ety2Colloc, etc are computed
 *          using analytical expressions for the tangent hyperbolic function.
 *
 *          Similarly, the staggered grids, \ref xStaggr, \ref yStaggr, and \ref zStaggr are initialized from the mid-point
 *          averaged values of the nodes in their corresponding transformed plane coordinates, \ref xi, \ref et, and \ref zt
 *          respectively using the tangent hyperbolic function.
 *          As before, the grid derivative terms for the staggered points are also computed using analytical expressions.
 *
 * \param   dim is an integer value that defines the direction along which tan-hyp grid is to be generated: 0 -> X, 1 -> Y, 2 -> Z
 ********************************************************************************************************************************************
 */
void grid::createTanHypGrid(int dim) {
    int i;

#ifndef TEST_RUN
    if (rankData.rank == 0) {
        switch (dim) {
            case 0: std::cout << "Generating tangent hyperbolic grid along X direction" << std::endl;
                    break;
            case 1: std::cout << "Generating tangent hyperbolic grid along Y direction" << std::endl;
                    break;
            case 2: std::cout << "Generating tangent hyperbolic grid along Z direction" << std::endl;
                    break;
        }
    }
#endif

    if (dim == 0) {
        // STAGGERED X-GRID POINTS FROM UNIFORM XI-GRID POINTS AND THEIR METRICS
        for (i = 0; i < coreSize(0); i++) {
            x(i) = xLen*(1.0 - tanh(thBeta[0]*(1.0 - 2.0*xi(i)))/tanh(thBeta[0]))/2.0;

            xi_x(i) = tanh(thBeta[0])/(thBeta[0]*xLen*(1.0 - pow((1.0 - 2.0*x(i)/xLen)*tanh(thBeta[0]), 2)));
            xixx(i) = -4.0*pow(tanh(thBeta[0]), 3)*(1.0 - 2.0*x(i)/xLen)/(thBeta[0]*xLen*xLen*pow(1.0 - pow(tanh(thBeta[0])*(1.0 - 2.0*x(i)/xLen), 2), 2));
            xix2(i) = pow(xi_x(i), 2.0);
        }
    }

#ifndef PLANAR
    if (dim == 1) {
        // STAGGERED Y-GRID POINTS FROM UNIFORM ETA-GRID POINTS AND THEIR METRICS
        for (i = 0; i < coreSize(1); i++) {
            y(i) = yLen*(1.0 - tanh(thBeta[1]*(1.0 - 2.0*et(i)))/tanh(thBeta[1]))/2.0;

            et_y(i) = tanh(thBeta[1])/(thBeta[1]*yLen*(1.0 - pow((1.0 - 2.0*y(i)/yLen)*tanh(thBeta[1]), 2)));
            etyy(i) = -4.0*pow(tanh(thBeta[1]), 3)*(1.0 - 2.0*y(i)/yLen)/(thBeta[1]*yLen*yLen*pow(1.0 - pow(tanh(thBeta[1])*(1.0 - 2.0*y(i)/yLen), 2), 2));
            ety2(i) = pow(et_y(i), 2.0);
        }
    }
#endif

    if (dim == 2) {
        // STAGGERED Z-GRID POINTS FROM UNIFORM ZETA-GRID POINTS AND THEIR METRICS
        for (i = 0; i < coreSize(2); i++) {
            z(i) = zLen*(1.0 - tanh(thBeta[2]*(1.0 - 2.0*zt(i)))/tanh(thBeta[2]))/2.0;

            zt_z(i) = tanh(thBeta[2])/(thBeta[2]*zLen*(1.0 - pow((1.0 - 2.0*z(i)/zLen)*tanh(thBeta[2]), 2)));
            ztzz(i) = -4.0*pow(tanh(thBeta[2]), 3)*(1.0 - 2.0*z(i)/zLen)/(thBeta[2]*zLen*zLen*pow(1.0 - pow(tanh(thBeta[2])*(1.0 - 2.0*z(i)/zLen), 2), 2));
            ztz2(i) = pow(zt_z(i), 2.0);
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronize grid data across the MPI sub-domains
 *
 *          The function updates the pad regions of the MPI sub-domains with grid data from adjacent sub-domains.
 *          Special MPI data structures are created temporarily for this one-off operation.
 *
 ********************************************************************************************************************************************
 */
void grid::syncGrid() {
    // Points at the left and right pads of the MPI-subdomains need to be taken from
    // neighbouring sub-domains, and this MPI datatype will be used for this purpose.
    MPI_Datatype padGrid;
    MPI_Status srStatus;
    int count;

    // At the leftmost and rightmost subdomains, the pad regions have to be populated differently.
    // The following blitz::Range objects identify the pad region and the corresponding regions
    // within the computational domain from which grid coordinates have to be pulled to appropriately
    // set the coordinates in the pads.
    blitz::Range xLftPad, xRgtPad, xLftPts, xRgtPts;
    blitz::Range yLftPad, yRgtPad, yLftPts, yRgtPts;
    blitz::Range zLftPad, zRgtPad, zLftPts, zRgtPts;

    // There is an assumption here that the pad widths are same along X and Y axes here
    // This is not a necessary assumption, but it allows me to write less code here.
    // Moreover, it is a pretty weird configuration to use different padwidths along X and Y axes
    count = padWidths(0);
    MPI_Type_contiguous(count, MPI_FP_REAL, &padGrid);
    MPI_Type_commit(&padGrid);

    // First, let us deal with the staggered grid
    // Along X Axis
    xLftPad = blitz::Range(-padWidths(0), -1, 1);
    xLftPts = blitz::Range(padWidths(0) , 1, -1);

    xRgtPad = blitz::Range(coreSize(0), coreSize(0) + padWidths(0) - 1, 1);
    xRgtPts = blitz::Range(coreSize(0) - 2, coreSize(0) - padWidths(0) - 1, -1);

    // Exchange the pad points across processors for staggered grid data along X axis
    MPI_Sendrecv(&x(1), 1, padGrid, rankData.faceRanks(0), 1, &x(coreSize(0)), 1, padGrid, rankData.faceRanks(1), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&x(coreSize(0) - padWidths(0) - 1), 1, padGrid, rankData.faceRanks(1), 2, &x(-padWidths(0)), 1, padGrid, rankData.faceRanks(0), 2, MPI_COMM_WORLD, &srStatus);

    MPI_Sendrecv(&xi_x(1), 1, padGrid, rankData.faceRanks(0), 1, &xi_x(coreSize(0)), 1, padGrid, rankData.faceRanks(1), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&xi_x(coreSize(0) - padWidths(0) - 1), 1, padGrid, rankData.faceRanks(1), 2, &xi_x(-padWidths(0)), 1, padGrid, rankData.faceRanks(0), 2, MPI_COMM_WORLD, &srStatus);

    MPI_Sendrecv(&xixx(1), 1, padGrid, rankData.faceRanks(0), 1, &xixx(coreSize(0)), 1, padGrid, rankData.faceRanks(1), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&xixx(coreSize(0) - padWidths(0) - 1), 1, padGrid, rankData.faceRanks(1), 2, &xixx(-padWidths(0)), 1, padGrid, rankData.faceRanks(0), 2, MPI_COMM_WORLD, &srStatus);

    MPI_Sendrecv(&xix2(1), 1, padGrid, rankData.faceRanks(0), 1, &xix2(coreSize(0)), 1, padGrid, rankData.faceRanks(1), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&xix2(coreSize(0) - padWidths(0) - 1), 1, padGrid, rankData.faceRanks(1), 2, &xix2(-padWidths(0)), 1, padGrid, rankData.faceRanks(0), 2, MPI_COMM_WORLD, &srStatus);

    // Whether domain is periodic or not, left-most point and right-most points will be differently calculated
    if (rankData.xRank == 0) {
        x(xLftPad) = -x(xLftPts);
        xixx(xLftPad) = -xixx(xLftPad);
    }
    if (rankData.xRank == rankData.npX - 1) {
        x(xRgtPad) = 2.0*xLen - x(xRgtPts);
        xixx(xRgtPad) = -xixx(xRgtPad);
    }

    // Along Y Axis
#ifndef PLANAR
    yLftPad = blitz::Range(-padWidths(1), -1, 1);
    yLftPts = blitz::Range(padWidths(1) , 1, -1);

    yRgtPad = blitz::Range(coreSize(1), coreSize(1) + padWidths(1) - 1, 1);
    yRgtPts = blitz::Range(coreSize(1) - 2, coreSize(1) - padWidths(1) - 1, -1);

    // Exchange the pad points across processors for staggered grid data along Y axis
    MPI_Sendrecv(&y(1), 1, padGrid, rankData.faceRanks(2), 1, &y(coreSize(1)), 1, padGrid, rankData.faceRanks(3), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&y(coreSize(1) - padWidths(1) - 1), 1, padGrid, rankData.faceRanks(3), 2, &y(-padWidths(1)), 1, padGrid, rankData.faceRanks(2), 2, MPI_COMM_WORLD, &srStatus);

    MPI_Sendrecv(&et_y(1), 1, padGrid, rankData.faceRanks(2), 1, &et_y(coreSize(1)), 1, padGrid, rankData.faceRanks(3), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&et_y(coreSize(1) - padWidths(1) - 1), 1, padGrid, rankData.faceRanks(3), 2, &et_y(-padWidths(1)), 1, padGrid, rankData.faceRanks(2), 2, MPI_COMM_WORLD, &srStatus);

    MPI_Sendrecv(&etyy(1), 1, padGrid, rankData.faceRanks(2), 1, &etyy(coreSize(1)), 1, padGrid, rankData.faceRanks(3), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&etyy(coreSize(1) - padWidths(1) - 1), 1, padGrid, rankData.faceRanks(3), 2, &etyy(-padWidths(1)), 1, padGrid, rankData.faceRanks(2), 2, MPI_COMM_WORLD, &srStatus);

    MPI_Sendrecv(&ety2(1), 1, padGrid, rankData.faceRanks(2), 1, &ety2(coreSize(1)), 1, padGrid, rankData.faceRanks(3), 1, MPI_COMM_WORLD, &srStatus);
    MPI_Sendrecv(&ety2(coreSize(1) - padWidths(1) - 1), 1, padGrid, rankData.faceRanks(3), 2, &ety2(-padWidths(1)), 1, padGrid, rankData.faceRanks(2), 2, MPI_COMM_WORLD, &srStatus);

    // Whether domain is periodic or not, front-most point and rear-most points will be differently calculated
    if (rankData.yRank == 0) {
        y(yLftPad) = -y(yLftPts);
        etyy(yLftPad) = -etyy(yLftPad);
    }
    if (rankData.yRank == rankData.npY - 1) {
        y(yRgtPad) = 2.0*yLen - y(yRgtPts);
        etyy(yRgtPad) = -etyy(yRgtPad);
    }
#endif

    // Along Z Axis
    zLftPad = blitz::Range(-padWidths(2), -1, 1);
    zLftPts = blitz::Range(padWidths(2) , 1, -1);

    zRgtPad = blitz::Range(coreSize(2), coreSize(2) + padWidths(2) - 1, 1);
    zRgtPts = blitz::Range(coreSize(2) - 2, coreSize(2) - padWidths(2) - 1, -1);

    // Whether domain is periodic or not, top-most point and bottom-most points will be differently calculated
    z(zLftPad) = -z(zLftPts);
    ztzz(zLftPad) = -ztzz(zLftPad);

    z(zRgtPad) = 2.0*zLen - z(zRgtPts);
    ztzz(zRgtPad) = -ztzz(zRgtPad);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to check the anisotropy of the grid
 *
 *          The function is called only if non-uniform grid is made.
 *          It scans through the entire grid cell-by-cell and checks the 2 aspect-ratios of the cell (one for 2D grids)
 *          The maximum value of aspect ratio is stored.
 *          Each MPI sub-domain checks within its limits and an MPI_Reduce call gets the global maximum.
 *
 ********************************************************************************************************************************************
 */
void grid::checkAnisotropy() {
    real cellMaxAR;
    real xWidth, zWidth;
    real localMax, globalMax;

#ifndef PLANAR
    real yWidth;
    real aRatio, bRatio;
    real xyRatio, yzRatio, zxRatio;
#endif

    localMax = 0.0;
#ifdef PLANAR
    for (int i = 0; i <= collocCoreDomain.ubound(0) + 1; i++) {
        for (int k = 0; k <= collocCoreDomain.ubound(2) + 1; k++) {
            xWidth = xColloc(i-1) - xColloc(i);
            zWidth = zColloc(k-1) - zColloc(k);
            cellMaxAR = std::max(xWidth/zWidth, zWidth/xWidth);
            if (cellMaxAR > localMax) localMax = cellMaxAR;

        }
    }
#else
    for (int i = 0; i <= collocCoreDomain.ubound(0) + 1; i++) {
        for (int j = 0; j <= collocCoreDomain.ubound(1) + 1; j++) {
            for (int k = 0; k <= collocCoreDomain.ubound(2) + 1; k++) {
                xWidth = xColloc(i-1) - xColloc(i);
                yWidth = yColloc(j-1) - yColloc(j);
                zWidth = zColloc(k-1) - zColloc(k);
                xyRatio = std::max(xWidth/yWidth, yWidth/xWidth);
                yzRatio = std::max(yWidth/zWidth, zWidth/yWidth);
                zxRatio = std::max(zWidth/xWidth, xWidth/zWidth);
                aRatio = std::max(xyRatio, yzRatio);
                bRatio = std::max(yzRatio, zxRatio);
                cellMaxAR = std::max(aRatio, bRatio);
                if (cellMaxAR > localMax) localMax = cellMaxAR;
            }
        }
    }
#endif

    MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (globalMax > 7.0) {
        if (rankData.rank == 0) std::cout << "\nWARNING: Grid anisotropy exceeds limits. Finite-difference calculations will be inaccurate" << std::endl;
    } else {
        if (rankData.rank == 0) std::cout << "\nMaximum grid anisotropy is " << globalMax << std::endl;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to gather global data about the grid into every rank
 *
 *          In certain cases, for example, implementation of non-homogeneous boundary conditions,
 *          the data regarding the global extents of the grid and the global coordinates are necessary.
 *          For this reason, the global grid data is made available to each rank through an MPI_Allgather call.
 *
 ********************************************************************************************************************************************
 */
void grid::gatherGlobal() {
    int i;
    int locSize, locDisp;
    int maxRank = std::max(rankData.npX, rankData.npY);

    int arrSize[maxRank];
    int arrDisp[maxRank];

    blitz::TinyVector<int, 3> globalSize;
    blitz::TinyVector<int, 3> globalReIndexVal;

    globalSize = globalSize + 2*padWidths;
    globalReIndexVal = -padWidths;

    xGlobal.resize(globalSize(0));     xGlobal.reindexSelf(globalReIndexVal(0));        xGlobal = 0.0;

#ifndef PLANAR
    yGlobal.resize(globalSize(1));     yGlobal.reindexSelf(globalReIndexVal(1));        yGlobal = 0.0;
#endif

    zGlobal.resize(globalSize(2));     zGlobal.reindexSelf(globalReIndexVal(2));        zGlobal = 0.0;

    // GATHERING THE STAGGERED GRID ALONG X-DIRECTION
    locSize = x.size() - 2*padWidths(0);
    locDisp = subarrayStarts(0);
    if (rankData.xRank == rankData.npX-1) {
        locSize += 2*padWidths(0);
    }
    MPI_Allgather(&locSize, 1, MPI_INT, arrSize, 1, MPI_INT, rankData.MPI_ROW_COMM);
    MPI_Allgather(&locDisp, 1, MPI_INT, arrDisp, 1, MPI_INT, rankData.MPI_ROW_COMM);
    MPI_Allgatherv(x.dataFirst(), locSize, MPI_FP_REAL, xGlobal.dataFirst(), arrSize, arrDisp, MPI_FP_REAL, rankData.MPI_ROW_COMM);

#ifndef PLANAR
    // GATHERING THE STAGGERED GRID ALONG Y-DIRECTION
    locSize = y.size() - 2*padWidths(1);
    locDisp = subarrayStarts(1);
    if (rankData.yRank == rankData.npY-1) {
        locSize += 2*padWidths(1);
    }
    MPI_Allgather(&locSize, 1, MPI_INT, arrSize, 1, MPI_INT, rankData.MPI_COL_COMM);
    MPI_Allgather(&locDisp, 1, MPI_INT, arrDisp, 1, MPI_INT, rankData.MPI_COL_COMM);
    MPI_Allgatherv(y.dataFirst(), locSize, MPI_FP_REAL, yGlobal.dataFirst(), arrSize, arrDisp, MPI_FP_REAL, rankData.MPI_COL_COMM);
#endif

    // GLOBAL AND LOCAL STAGGERED GRIDS ALONG Z-DIRECTION ARE SAME FOR ALL RANKS
    for (i = -padWidths(2); i < globalSize(2) + padWidths(2); i++) {
        zGlobal(i) = z(i);
    }
}
