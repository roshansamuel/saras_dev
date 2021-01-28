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
/*! \file hydro.cc
 *
 *  \brief Definitions of common functions for both 2D and 3D runs of the solver class hydro - this class solves the basic Navier-Stokes equation.
 *  \sa hydro.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "hydro.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the base hydro class
 *
 *          The short base constructor of the hydro class merely assigns the const references to the grid and parser
 *          class instances being used in the solver.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
hydro::hydro(const grid &mesh, const parser &solParam, parallel &mpiParam):
            V(mesh, "V"),
            P(mesh, "P"),
            mesh(mesh),
            inputParams(solParam),
            mpiData(mpiParam) { }


/**
 ********************************************************************************************************************************************
 * \brief   The core publicly accessible function of the \ref hydro class to solve the Navier-Stokes equations
 *
 *          The NSE are integrated in time from within this function by calling \ref hydro#timeAdvance in a loop.
 *          The function keeps track of the non-dimensional time with \ref time and number of iterations with \ref iterCount.
 *          Both these values are continuously incremented from within the loop, and finally, when \ref time has reached the
 *          user-ser value in \ref parser#tMax "tMax", the time-integration loop is broken and the program exits.
 ********************************************************************************************************************************************
 */
void hydro::solvePDE() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to enable/disable periodic data transfer as per the problem
 *
 *          The function checks the xPer, yPer and zPer flags in the parser class
 *          and enables/disables MPI data transfer at boundaries accordingly
 *          By default, the MPI neighbours at boundaries are set for periodic data-transfer.
 *          This has to be disabled if the problem has non-periodic boundaries.
 *
 ********************************************************************************************************************************************
 */
void hydro::checkPeriodic() {
    // Disable periodic data transfer by setting neighbouring ranks of boundary sub-domains to NULL
    // Left and right walls
    if (not inputParams.xPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along X Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.xRank == 0) {
            mpiData.faceRanks(0) = MPI_PROC_NULL;

            mpiData.edgeRanks(0) = MPI_PROC_NULL;
            mpiData.edgeRanks(1) = MPI_PROC_NULL;
        }

        if (mpiData.xRank == mpiData.npX-1) {
            mpiData.faceRanks(1) = MPI_PROC_NULL;

            mpiData.edgeRanks(2) = MPI_PROC_NULL;
            mpiData.edgeRanks(3) = MPI_PROC_NULL;
        }
    }

    // Front and rear walls
#ifdef PLANAR
    // Front and rear walls are by default non-periodic for 2D simulations
    if (mpiData.yRank == 0)             mpiData.faceRanks(2) = MPI_PROC_NULL;
    if (mpiData.yRank == mpiData.npY-1) mpiData.faceRanks(3) = MPI_PROC_NULL;

#else
    if (not inputParams.yPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Y Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.yRank == 0) {
            mpiData.faceRanks(2) = MPI_PROC_NULL;

            mpiData.edgeRanks(0) = MPI_PROC_NULL;
            mpiData.edgeRanks(2) = MPI_PROC_NULL;
        }

        if (mpiData.yRank == mpiData.npY-1) {
            mpiData.faceRanks(3) = MPI_PROC_NULL;

            mpiData.edgeRanks(1) = MPI_PROC_NULL;
            mpiData.edgeRanks(3) = MPI_PROC_NULL;
        }
    }
#endif

    // Inform user about BC along top and bottom walls
    if (not inputParams.zPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Z Direction" << std::endl;
            std::cout << std::endl;
        }
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the forcing terms for velocity
 *
 *          The forcing terms for the velocity field are initialized here.
 *          Out of the different forcings available in the force class,
 *          the appropriate forcing is chosen according to the parameters set by the user.
 ********************************************************************************************************************************************
 */
void hydro::initVForcing() {
    switch (inputParams.forceType) {
        case 0:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with zero velocity forcing" << std::endl << std::endl;
            V.vForcing = new zeroForcing(mesh, V);
            break;
        case 1:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with random velocity forcing" << std::endl << std::endl;
            V.vForcing = new randomForcing(mesh, V);
            break;
        case 2:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with rotation" << std::endl << std::endl;
            V.vForcing = new coriolisForce(mesh, V);
            break;
        case 5:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with constant pressure gradient along X-direction" << std::endl << std::endl;
            V.vForcing = new constantPGrad(mesh, V);
            break;
        default:
            if (mpiData.rank == 0) std::cout << "WARNING: Chosen velocity forcing is incompatible with hydrodynamics runs. Defaulting to zero forcing" << std::endl << std::endl;
            V.vForcing = new zeroForcing(mesh, V);
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the boundary conditions for velocity
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::initVBCs() {
    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        V.uLft = new dirichletFC(mesh, V.Vx, 0, 1.0);
        V.uRgt = new neumannFC(mesh, V.Vx, 1, 0.0);
    } else {
        // NO-PENETRATION BCS
        V.uLft = new dirichletFC(mesh, V.Vx, 0, 0.0);
        V.uRgt = new dirichletFC(mesh, V.Vx, 1, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    V.uFrn = new dirichletCC(mesh, V.Vx, 2, 0.0);
    V.uBak = new dirichletCC(mesh, V.Vx, 3, 0.0);
#endif

    if (inputParams.zPer) {
        // PERIODIC BC
        V.uBot = new periodicCC(mesh, V.Vx, 4);
        V.uTop = new periodicCC(mesh, V.Vx, 5);
    } else {
        if (inputParams.probType == 1) {
            // NO-SLIP BCS FOR LDC
            V.uBot = new dirichletCC(mesh, V.Vx, 4, 0.0);
            V.uTop = new dirichletCC(mesh, V.Vx, 5, 1.0);
        } else {
            // NO-SLIP BCS
            V.uBot = new dirichletCC(mesh, V.Vx, 4, 0.0);
            V.uTop = new dirichletCC(mesh, V.Vx, 5, 0.0);
        }
    }

#ifndef PLANAR
    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        V.vLft = new dirichletCC(mesh, V.Vy, 0, 0.0);
        V.vRgt = new neumannCC(mesh, V.Vy, 1, 0.0);
    } else {
        // NO-SLIP BCS
        V.vLft = new dirichletCC(mesh, V.Vy, 0, 0.0);
        V.vRgt = new dirichletCC(mesh, V.Vy, 1, 0.0);
    }

    // NO-PENETRATION BCS
    V.vFrn = new dirichletFC(mesh, V.Vy, 2, 0.0);
    V.vBak = new dirichletFC(mesh, V.Vy, 3, 0.0);

    if (inputParams.zPer) {
        // PERIODIC BC
        V.vBot = new periodicCC(mesh, V.Vy, 4);
        V.vTop = new periodicCC(mesh, V.Vy, 5);
    } else {
        // NO-SLIP BCS
        V.vBot = new dirichletCC(mesh, V.Vy, 4, 0.0);
        V.vTop = new dirichletCC(mesh, V.Vy, 5, 0.0);
    }
#endif

    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        V.wLft = new dirichletCC(mesh, V.Vz, 0, 0.0);
        V.wRgt = new neumannCC(mesh, V.Vz, 1, 0.0);
    } else {
        // NO-SLIP BCS
        V.wLft = new dirichletCC(mesh, V.Vz, 0, 0.0);
        V.wRgt = new dirichletCC(mesh, V.Vz, 1, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    V.wFrn = new dirichletCC(mesh, V.Vz, 2, 0.0);
    V.wBak = new dirichletCC(mesh, V.Vz, 3, 0.0);
#endif

    if (inputParams.zPer) {
        // PERIODIC BC
        V.wBot = new periodicFC(mesh, V.Vz, 4);
        V.wTop = new periodicFC(mesh, V.Vz, 5);
    } else {
        // NO-SLIP BCS
        V.wBot = new dirichletFC(mesh, V.Vz, 4, 0.0);
        V.wTop = new dirichletFC(mesh, V.Vz, 5, 0.0);
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the boundary conditions for pressure
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::initPBCs() {
    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        P.tLft = new nullBC(mesh, P.F, 0);
        P.tRgt = new neumannCC(mesh, P.F, 1, 0.0);
    } else {
        // NO BCS REQUIRED FOR P WHEN V IS SPECIFED
        P.tLft = new nullBC(mesh, P.F, 0);
        P.tRgt = new nullBC(mesh, P.F, 1);
    }

#ifndef PLANAR
    // NO BCS REQUIRED FOR P WHEN V IS SPECIFED ON NO-SLIP WALLS
    P.tFrn = new nullBC(mesh, P.F, 2);
    P.tBak = new nullBC(mesh, P.F, 3);
#endif

    if (inputParams.zPer) {
        // PERIODIC BC
        P.tBot = new periodicCC(mesh, P.F, 4);
        P.tTop = new periodicCC(mesh, P.F, 5);
    } else {
        // NO BCS REQUIRED FOR P WHEN V IS SPECIFED
        P.tBot = new nullBC(mesh, P.F, 4);
        P.tTop = new nullBC(mesh, P.F, 5);
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether periodic BC is being implemented properly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls imposeUBCs, imposeVBCs and imposeWBCs and checks if the correct values of the functions are imposed at boundaries
 ********************************************************************************************************************************************
 */
real hydro::testPeriodic() { return 0; };
