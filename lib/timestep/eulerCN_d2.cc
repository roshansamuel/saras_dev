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
/*! \file eulerCN_d2.cc
 *
 *  \brief Definitions for functions of class timestep
 *  \sa timestep.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "timestep.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the timestep class
 *
 *          The empty constructer merely initializes the local reference to the global mesh variable.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 ********************************************************************************************************************************************
 */
eulerCN_d2::eulerCN_d2(const grid &mesh, const real &dt, vfield &V, sfield &P):
    timestep(mesh, dt, V, P),
    mgSolver(mesh, mesh.inputParams)
{
    setCoefficients();

    maxIterations = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded constructor of the timestep class
 *
 *          In addition to the variables for hydrodynamics simulations, variables for scalar simulations are also intialized
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 ********************************************************************************************************************************************
 */
eulerCN_d2::eulerCN_d2(const grid &mesh, const real &dt, vfield &V, sfield &P, sfield &T):
    timestep(mesh, dt, V, P, T),
    mgSolver(mesh, mesh.inputParams)
{
    setCoefficients();

    maxIterations = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the solution using Euler method and Implicit Crank-Nicholson method
 *
 *          The non-linear terms are advanced using explicit Euler method, while the duffusion terms are
 *          advanced by semi-implicit Crank-Nicholson method.
 *
 ********************************************************************************************************************************************
 */
void eulerCN_d2::timeAdvance(vfield &V, sfield &P) {
    static plainvf nseRHS(mesh, V);

    nseRHS = 0.0;

    // First compute the explicit part of the semi-implicit viscous term and divide it by Re
    V.computeDiff(nseRHS);
    nseRHS *= nu;

    // Compute the non-linear term and subtract it from the RHS
    V.computeNLin(V, nseRHS);

    // Add the velocity forcing term
    V.vForcing->addForcing(nseRHS);

    // Subtract the pressure gradient term
    pressureGradient = 0.0;
    P.gradient(pressureGradient, V);
    nseRHS -= pressureGradient;

    // Multiply the entire RHS with dt and add the velocity of previous time-step to advance by explicit Euler method
    nseRHS *= dt;
    nseRHS += V;

    // Synchronize the RHS term across all processors by updating its sub-domain pads
    nseRHS.syncData();

    // Using the RHS term computed, compute the guessed velocity of CN method iteratively (and store it in V)
    solveVx(V, nseRHS);
    solveVz(V, nseRHS);

    // Calculate the rhs for the poisson solver (mgRHS) using the divergence of guessed velocity in V
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;

    // IF THE POISSON SOLVER IS BEING TESTED, THE RHS IS SET TO ONE.
    // THIS IS FOR TESTING ONLY AND A SINGLE TIME ADVANCE IS PERFORMED IN THIS TEST
#ifdef TEST_POISSON
    mgRHS.F = 1.0;
#endif

    // Using the calculated mgRHS, evaluate pressure correction (Pp) using multi-grid method
    mgSolver.mgSolve(Pp, mgRHS);

    // Synchronise the pressure correction term across processors
    Pp.syncData();

    // IF THE POISSON SOLVER IS BEING TESTED, THE PRESSURE IS SET TO ZERO.
    // THIS WAY, AFTER THE SOLUTION OF MG SOLVER, Pp, IS DIRECTLY WRITTEN INTO P AND AVAILABLE FOR PLOTTING
    // THIS IS FOR TESTING ONLY AND A SINGLE TIME ADVANCE IS PERFORMED IN THIS TEST
#ifdef TEST_POISSON
    P.F = 0.0;
#endif

    // Add the pressure correction term to the pressure field of previous time-step, P
    P += Pp;

    // Finally get the velocity field at end of time-step by subtracting the gradient of pressure correction from V
    Pp.gradient(pressureGradient, V);
    pressureGradient *= dt;
    V -= pressureGradient;

    // Impose boundary conditions on the updated velocity field, V
    V.imposeBCs();
}


void eulerCN_d2::solveVx(vfield &V, plainvf &nseRHS) {
    int iterCount = 0;
    real maxError = 0.0;
    static blitz::Array<real, 3> tempVx(V.Vx.F.lbound(), V.Vx.F.shape());

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(iY) shared(V) shared(nseRHS) shared(tempVx)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                tempVx(iX, iY, iZ) = ((hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                       hx2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                        dt * nu / ( hz2hx2 * 2.0) + nseRHS.Vx(iX, iY, iZ)) /
                 (1.0 + dt * nu * ((hz2 * mesh.xix2Colloc(iX) + hx2 * mesh.ztz2Staggr(iZ)))/hz2hx2);
            }
        }

        V.Vx.F = tempVx;

        V.imposeVxBC();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(iY) shared(V) shared(tempVx)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                tempVx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - 0.5 * dt * nu * (
                          mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx2) +
                          mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz2));
            }
        }

        tempVx(V.Vx.fBulk) = abs(tempVx(V.Vx.fBulk) - nseRHS.Vx(V.Vx.fBulk));

        maxError = blitz::max(tempVx);

        if (maxError < mesh.inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vx not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void eulerCN_d2::solveVz(vfield &V, plainvf &nseRHS) {
    int iterCount = 0;
    real maxError = 0.0;
    static blitz::Array<real, 3> tempVz(V.Vz.F.lbound(), V.Vz.F.shape());

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(iY) shared(V) shared(nseRHS) shared(tempVz)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                tempVz(iX, iY, iZ) = ((hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                       hx2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                        dt * nu / ( hz2hx2 * 2.0) + nseRHS.Vz(iX, iY, iZ)) /
                 (1.0 + dt * nu * ((hz2 * mesh.xix2Staggr(iX) + hx2 * mesh.ztz2Colloc(iZ)))/hz2hx2);
            }
        }

        V.Vz.F = tempVz;

        V.imposeVzBC();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(iY) shared(V) shared(tempVz)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                tempVz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - 0.5 * dt * nu * (
                          mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx2) +
                          mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz2));
            }
        }

        tempVz(V.Vz.fBulk) = abs(tempVz(V.Vz.fBulk) - nseRHS.Vz(V.Vz.fBulk));

        maxError = blitz::max(tempVz);

        if (maxError < mesh.inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vz not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void eulerCN_d2::solveT(sfield &T, plainsf &tmpRHS) {
    int iterCount = 0;
    real maxError = 0.0;
    static blitz::Array<real, 3> tempT(T.F.F.lbound(), T.F.F.shape());

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(iY) shared(T) shared(tmpRHS) shared(tempT)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                tempT(iX, iY, iZ) = ((hz2 * mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) + T.F.F(iX-1, iY, iZ)) +
                                      hx2 * mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) + T.F.F(iX, iY, iZ-1))) *
                    dt * kappa / ( hz2hx2 * 2.0) + tmpRHS.F(iX, iY, iZ)) /
             (1.0 + dt * kappa * ((hz2 * mesh.xix2Staggr(iX) + hx2 * mesh.ztz2Colloc(iZ)))/hz2hx2);
            }
        }

        T.F.F = tempT;

        T.imposeBCs();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(iY) shared(T) shared(tempT)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                tempT(iX, iY, iZ) = T.F.F(iX, iY, iZ) - 0.5 * dt * kappa * (
                       mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX-1, iY, iZ)) / (hx2) +
                       mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY, iZ-1)) / (hz2));
            }
        }

        tempT(T.F.fBulk) = abs(tempT(T.F.fBulk) - tmpRHS.F(T.F.fBulk));

        maxError = blitz::max(tempT);

        if (maxError < mesh.inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of T not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the coefficients used for solving the implicit equations of U, V and W
 *
 *          The function assigns values to the variables \ref hx, \ref hy, etc.
 *          These coefficients are repeatedly used at many places in the iterative solver for implicit calculation of velocities.
 ********************************************************************************************************************************************
 */
void eulerCN_d2::setCoefficients() {
    hx2 = pow(mesh.dXi, 2.0);
    hz2 = pow(mesh.dZt, 2.0);

    hz2hx2 = pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);
};