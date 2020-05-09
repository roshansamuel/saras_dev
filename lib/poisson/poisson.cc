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
}


/**
 ********************************************************************************************************************************************
 * \brief   The core, publicly accessible function of poisson to compute the solution for the Poisson equation
 *
 ********************************************************************************************************************************************
 */
void poisson::mgSolve(plainsf &inFn, const plainsf &rhs) {
    for (int i=0; i <= inputParams.vcDepth; i++) {
        pressureData(i) = 0.0;
        residualData(i) = 0.0;

        // This was originally inside the V-Cycle loop. Check if this line being here is okay.
        smoothedPres(i) = 0.0;
    }

    // TRANSFER DATA FROM THE INPUT SCALAR FIELDS INTO THE DATA-STRUCTURES USED BY poisson
    residualData(0)(stagCore(0)) = rhs.F(stagCore(0));
    pressureData(0)(stagCore(0)) = inFn.F(stagCore(0));

    // PERFORM V-CYCLES AS MANY TIMES AS REQUIRED
    for (int i=0; i<inputParams.vcCount; i++) {
        vCycle();

        if (inputParams.mgError) {
            real mgResidual = computeError(inputParams.mgError);
            if (mesh.rankData.rank == 0) std::cout << std::endl << "Residual after V Cycle " << i << " is " << mgResidual << std::endl;
        }
    }

    // RETURN CALCULATED PRESSURE DATA
    inFn.F = pressureData(0)(blitz::RectDomain<3>(inFn.F.lbound(), inFn.F.ubound()));

#ifdef TEST_POISSON
    if (mesh.rankData.rank == 0) {
        real xDist, zDist;
        blitz::Array<real, 3> pAnalytic, tempArray;

        pAnalytic.resize(blitz::TinyVector<int, 3>(stagCore(0).ubound(0) - stagCore(0).lbound(0) + 1,
                                                   stagCore(0).ubound(1) - stagCore(0).lbound(1) + 1,
                                                   stagCore(0).ubound(2) - stagCore(0).lbound(2) + 1));
        pAnalytic.reindexSelf(blitz::TinyVector<int, 3>(stagCore(0).lbound(0),
                                                        stagCore(0).lbound(1),
                                                        stagCore(0).lbound(2)));
        pAnalytic = 0.0;

        for (int i=stagCore(0).lbound(0); i<=stagCore(0).ubound(0); i++) {
            xDist = hx(0)*(i - stagCore(0).ubound(0)/2);
            for (int k=stagCore(0).lbound(2); k<=stagCore(0).ubound(2); k++) {
                zDist = hz(0)*(k - stagCore(0).ubound(2)/2);

                pAnalytic(i, 0, k) = (xDist*xDist + zDist*zDist)/4.0;
            }
        }

        tempArray.resize(pAnalytic.shape());
        tempArray.reindexSelf(pAnalytic.lbound());

        tempArray = pAnalytic - pressureData(0)(stagCore(0));

        std::cout << std::endl;
        std::cout << "Maximum absolute deviation from analytic solution is: " << blitz::max(fabs(tempArray)) << std::endl;
        std::cout << std::endl;
    }
#endif
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to perform one loop of V-cycle
 *
 ********************************************************************************************************************************************
 */
void poisson::vCycle() {
    vLevel = 0;
    zeroBC = false;

    smooth(inputParams.preSmooth);

    zeroBC = true;

    for (int i=0; i<inputParams.vcDepth; i++) {
        computeResidual();
        smoothedPres(vLevel) = pressureData(vLevel);

        coarsen();
        pressureData(vLevel) = 0.0;

        (vLevel == inputParams.vcDepth)? solve(): smooth(inputParams.preSmooth);
    }

    for (int i=0; i<inputParams.vcDepth; i++) {
        prolong();

        pressureData(vLevel) += smoothedPres(vLevel);

        (vLevel == 0)? zeroBC = false: zeroBC = true;
        smooth(inputParams.postSmooth);
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the arrays used in multi-grid
 *
 ********************************************************************************************************************************************
 */
void poisson::initializeArrays() {
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
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the RectDomain variables for all future references throughout the poisson solver
 *
 ********************************************************************************************************************************************
 */
void poisson::setStagBounds() {
    blitz::TinyVector<int, 3> loBound, upBound;

    stagFull.resize(inputParams.vcDepth + 1);
    stagCore.resize(inputParams.vcDepth + 1);

    xEnd.resize(inputParams.vcDepth + 1);
    yEnd.resize(inputParams.vcDepth + 1);
    zEnd.resize(inputParams.vcDepth + 1);

    for (int i=0; i<=inputParams.vcDepth; i++) {
        // LOWER BOUND AND UPPER BOUND OF STAGGERED CORE - USED TO CONSTRUCT THE CORE SLICE
        loBound = 0, 0, 0;
#ifdef PLANAR
        upBound = mgSizeArray(localSizeIndex(0) - i) - 1, 0, mgSizeArray(localSizeIndex(2) - i) - 1;
#else
        upBound = mgSizeArray(localSizeIndex(0) - i) - 1, mgSizeArray(localSizeIndex(1) - i) - 1, mgSizeArray(localSizeIndex(2) - i) - 1;
#endif
        stagCore(i) = blitz::RectDomain<3>(loBound, upBound);

        // LOWER BOUND AND UPPER BOUND OF STAGGERED FULL SUB-DOMAIN - USED TO CONSTRUCT THE FULL SUB-DOMAIN SLICE
        loBound = -1, -1, -1;
        upBound = stagCore(i).ubound() - loBound;
        stagFull(i) = blitz::RectDomain<3>(loBound, upBound);

        // SET THE LIMTS FOR ARRAY LOOPS IN smooth FUNCTION, AND A FEW OTHER PLACES
        // WARNING: THESE VARIABLES HAVE SO FAR BEEN IMPLEMENTED ONLY IN smooth AND vCycle.
        // THE TEST FUNCTIONS HAVE NOT YET BEEN UPDATED WITH THESE
        xEnd(i) = stagCore(i).ubound(0);
#ifndef PLANAR
        yEnd(i) = stagCore(i).ubound(1);
#endif
        zEnd(i) = stagCore(i).ubound(2);
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the local size indices of sub-domains after MPI domain decomposition
 *
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
 ********************************************************************************************************************************************
 */
void poisson::copyStaggrDerivs() {
    xixx.resize(inputParams.vcDepth + 1);
    xix2.resize(inputParams.vcDepth + 1);
#ifndef PLANAR
    etyy.resize(inputParams.vcDepth + 1);
    ety2.resize(inputParams.vcDepth + 1);
#endif
    ztzz.resize(inputParams.vcDepth + 1);
    ztz2.resize(inputParams.vcDepth + 1);

    // WARNING:: BELOW METHOD IS WRONG! CORRECT THE STRIDES FOR EACH RANGE
    for(int i=0; i<=inputParams.vcDepth; i++) {
        xixx(i).resize(stagFull(i).ubound(0) - stagFull(i).lbound(0) + 1);
        xixx(i).reindexSelf(stagFull(i).lbound(0));
        xixx(i) = 0.0;
        xixx(i)(blitz::Range(0, stagCore(i).ubound(0), 1)) = mesh.xixxStaggr(blitz::Range(0, stagCore(i).ubound(0), 1));

        xix2(i).resize(stagFull(i).ubound(0) - stagFull(i).lbound(0) + 1);
        xix2(i).reindexSelf(stagFull(i).lbound(0));
        xix2(i) = 0.0;
        xix2(i)(blitz::Range(0, stagCore(i).ubound(0), 1)) = mesh.xix2Staggr(blitz::Range(0, stagCore(i).ubound(0), 1));

#ifndef PLANAR
        etyy(i).resize(stagFull(i).ubound(1) - stagFull(i).lbound(1) + 1);
        etyy(i).reindexSelf(stagFull(i).lbound(1));
        etyy(i) = 0.0;
        etyy(i)(blitz::Range(0, stagCore(i).ubound(1), 1)) = mesh.etyyStaggr(blitz::Range(0, stagCore(i).ubound(1), 1));

        ety2(i).resize(stagFull(i).ubound(1) - stagFull(i).lbound(1) + 1);
        ety2(i).reindexSelf(stagFull(i).lbound(1));
        ety2(i) = 0.0;
        ety2(i)(blitz::Range(0, stagCore(i).ubound(1), 1)) = mesh.ety2Staggr(blitz::Range(0, stagCore(i).ubound(1), 1));
#endif

        ztzz(i).resize(stagFull(i).ubound(2) - stagFull(i).lbound(2) + 1);
        ztzz(i).reindexSelf(stagFull(i).lbound(2));
        ztzz(i) = 0.0;
        ztzz(i)(blitz::Range(0, stagCore(i).ubound(2), 1)) = mesh.ztzzStaggr(blitz::Range(0, stagCore(i).ubound(2), 1));

        ztz2(i).resize(stagFull(i).ubound(2) - stagFull(i).lbound(2) + 1);
        ztz2(i).reindexSelf(stagFull(i).lbound(2));
        ztz2(i) = 0.0;
        ztz2(i)(blitz::Range(0, stagCore(i).ubound(2), 1)) = mesh.ztz2Staggr(blitz::Range(0, stagCore(i).ubound(2), 1));
    }
};


void poisson::coarsen() { };


void poisson::prolong() { };


void poisson::computeResidual() { };


void poisson::smooth(const int smoothCount) { };


real poisson::computeError(const int normOrder) {  return 0.0; };


void poisson::imposeBC() { };


poisson::~poisson() {
#ifdef TIME_RUN
    if (mesh.rankData.rank == 0) {
        std::cout << std::left << std::setw(50) << "Time taken in computation within smooth: "           << std::fixed << std::setprecision(6) << smothTimeComp << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in data-transfer within smooth: "         << std::fixed << std::setprecision(6) << smothTimeTran << std::endl;
    }
#endif
};
