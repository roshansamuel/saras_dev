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
/*! \file eulerCN_d3.cc
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
 *          The empty constructor merely initializes the local reference to the global mesh variable.
 *          Also, the maximum allowable number of iterations for the Jacobi iterative solver being used to solve for the
 *          velocities implicitly is set as \f$ N_{max} = N_x \times N_y \times N_z \f$, where \f$N_x\f$, \f$N_y\f$ and \f$N_z\f$
 *          are the number of grid points in the collocated grid at the local sub-domains along x, y and z directions
 *          respectively.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 ********************************************************************************************************************************************
 */
eulerCN_d3::eulerCN_d3(const grid &mesh, const real &dt, vfield &V, sfield &P):
    timestep(mesh, dt, V, P),
    mgSolver(mesh, mesh.inputParams)
{
    setCoefficients();

    // This upper limit on max iterations is an arbitrarily chosen function.
    // Using Nx x Ny x Nz as the upper limit may cause the run to freeze for very long time.
    // This can eat away a lot of core hours unnecessarily.
    // It remains to be seen if this upper limit is safe.
    maxIterations = int(std::pow(std::log(mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2)), 3));
    //maxIterations = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);

    // If LES switch is enabled, initialize LES model
    if (mesh.inputParams.lesModel) {
        if (mesh.rankData.rank == 0) {
            std::cout << "LES Switch is ON. Using stretched spiral vortex LES Model\n" << std::endl;
        }
        sgsLES = new spiral(mesh);

        Txx = new sfield(mesh, "Txx");
        Tyy = new sfield(mesh, "Tyy");
        Tzz = new sfield(mesh, "Tzz");
        Txy = new sfield(mesh, "Txy");
        Tyz = new sfield(mesh, "Tyz");
        Tzx = new sfield(mesh, "Tzx");

        qX = new sfield(mesh, "qX");
        qY = new sfield(mesh, "qY");
        qZ = new sfield(mesh, "qZ");

        Vxcc = new sfield(mesh, "Ucc");
        Vycc = new sfield(mesh, "Vcc");
        Vzcc = new sfield(mesh, "Wcc");
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to advance the solution using Euler method and Implicit Crank-Nicholson method
 *
 *          The non-linear terms are advanced using explicit Euler method, while the duffusion terms are
 *          advanced by semi-implicit Crank-Nicholson method.
 *          This overloaded function advances velocity and pressure fields for hydrodynamics simulations.
 *
 ********************************************************************************************************************************************
 */
void eulerCN_d3::timeAdvance(vfield &V, sfield &P) {
    static plainvf nseRHS(mesh, V);

    nseRHS = 0.0;

    // First compute the explicit part of the semi-implicit viscous term and divide it by Re
    V.computeDiff(nseRHS);
    nseRHS *= nu;

    // Compute the non-linear term and subtract it from the RHS
    V.computeNLin(V, nseRHS);

    // Add the velocity forcing term
    V.vForcing->addForcing(nseRHS);

    // Add sub-grid stress contribution from LES Model, if enabled
    if (mesh.inputParams.lesModel) {
        // Array limits for loops
        int xS, xE, yS, yE, zS, zE;

        // Temporary variables to store output from spiral LES solver
        double sTxx, sTyy, sTzz, sTxy, sTyz, sTzx;

        // These are three 3x3x3 arrays containing local interpolated velocities
        // These are used to calculate the structure function within the spiral les routine
        blitz::Array<double, 3> u(3, 3, 3), v(3, 3, 3), w(3, 3, 3);

        // The following 3x3 matrix stores the velocity gradient tensor
        blitz::Array<double, 2> dudx(3, 3);

        // These 9 arrays store components of the velocity gradient tensor intially
        // Then they are reused to store the derivatives of stress tensor to calculate its divergence
        blitz::Array<double, 3> A11, A12, A13;
        blitz::Array<double, 3> A21, A22, A23;
        blitz::Array<double, 3> A31, A32, A33;

        // The 9 arrays of tensor components have the same dimensions and limits as a cell centered variable
        A11.resize(P.F.fSize);      A11.reindexSelf(P.F.flBound);
        A12.resize(P.F.fSize);      A12.reindexSelf(P.F.flBound);
        A13.resize(P.F.fSize);      A13.reindexSelf(P.F.flBound);
        A21.resize(P.F.fSize);      A21.reindexSelf(P.F.flBound);
        A22.resize(P.F.fSize);      A22.reindexSelf(P.F.flBound);
        A23.resize(P.F.fSize);      A23.reindexSelf(P.F.flBound);
        A31.resize(P.F.fSize);      A31.reindexSelf(P.F.flBound);
        A32.resize(P.F.fSize);      A32.reindexSelf(P.F.flBound);
        A33.resize(P.F.fSize);      A33.reindexSelf(P.F.flBound);

        // First interpolate all velocities to cell centers
        // Then U, V and W data are available at all cell centers except ghost points
        P.interTempF = 0.0;
        for (unsigned int i=0; i < P.F.VxIntSlices.size(); i++) {
            P.interTempF(P.F.fCore) += V.Vx.F(P.F.VxIntSlices(i));
        }
        Vxcc->F.F = P.interTempF/P.F.VxIntSlices.size();

        P.interTempF = 0.0;
        for (unsigned int i=0; i < P.F.VyIntSlices.size(); i++) {
            P.interTempF(P.F.fCore) += V.Vy.F(P.F.VyIntSlices(i));
        }
        Vycc->F.F = P.interTempF/P.F.VyIntSlices.size();

        P.interTempF = 0.0;
        for (unsigned int i=0; i < P.F.VzIntSlices.size(); i++) {
            P.interTempF(P.F.fCore) += V.Vz.F(P.F.VzIntSlices(i));
        }
        Vzcc->F.F = P.interTempF/P.F.VzIntSlices.size();

        Vxcc->derS.calcDerivative1_x(A11);
        Vxcc->derS.calcDerivative1_y(A12);
        Vxcc->derS.calcDerivative1_z(A13);
        Vycc->derS.calcDerivative1_x(A21);
        Vycc->derS.calcDerivative1_y(A22);
        Vycc->derS.calcDerivative1_z(A23);
        Vzcc->derS.calcDerivative1_x(A31);
        Vzcc->derS.calcDerivative1_y(A32);
        Vzcc->derS.calcDerivative1_z(A33);

        // Use only cell centers like a collocated grid and compute T tensor
        // Since U, V, and W data is available only in the core,
        // Adjust the loop limits so that the boundary points are excluded
        // to compute derivatives and structure functions correctly.
        xS = P.F.fCore.lbound(0) + 1; xE = P.F.fCore.ubound(0) - 1;
        yS = P.F.fCore.lbound(1) + 1; yE = P.F.fCore.ubound(1) - 1;
        zS = P.F.fCore.lbound(2) + 1; zE = P.F.fCore.ubound(2) - 1;

        for (int iX = xS; iX <= xE; iX++) {
            double dx = mesh.xColloc(iX - 1) - mesh.xColloc(iX);
            for (int iY = yS; iY <= yE; iY++) {
                double dy = mesh.yColloc(iY - 1) - mesh.yColloc(iY);
                for (int iZ = zS; iZ <= zE; iZ++) {
                    double dz = mesh.zColloc(iZ - 1) - mesh.zColloc(iZ);

                    // Cutoff wavelength
                    double del = std::cbrt(dx*dy*dz);

                    u = Vxcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                    v = Vycc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                    w = Vzcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));

                    // A 3 x 3 matrix of velocity gradient
                    dudx = A11(iX, iY, iZ), A12(iX, iY, iZ), A13(iX, iY, iZ),
                           A21(iX, iY, iZ), A22(iX, iY, iZ), A23(iX, iY, iZ),
                           A31(iX, iY, iZ), A32(iX, iY, iZ), A33(iX, iY, iZ);

                    double x[3] = {mesh.xStaggr(iX-1), mesh.xStaggr(iX), mesh.xStaggr(iX+1)};
                    double y[3] = {mesh.yStaggr(iY-1), mesh.yStaggr(iY), mesh.yStaggr(iY+1)};
                    double z[3] = {mesh.zStaggr(iZ-1), mesh.zStaggr(iZ), mesh.zStaggr(iZ+1)};

                    sgsLES->sgs_stress(u, v, w, dudx, x, y, z, nu, del,
                                    &sTxx, &sTyy, &sTzz, &sTxy, &sTyz, &sTzx);

                    Txx->F.F(iX, iY, iZ) = sTxx;
                    Tyy->F.F(iX, iY, iZ) = sTyy;
                    Tzz->F.F(iX, iY, iZ) = sTzz;
                    Txy->F.F(iX, iY, iZ) = sTxy;
                    Tyz->F.F(iX, iY, iZ) = sTyz;
                    Tzx->F.F(iX, iY, iZ) = sTzx;
                }
            }
        }

        Txx->syncData();
        Tyy->syncData();
        Tzz->syncData();
        Txy->syncData();
        Tyz->syncData();
        Tzx->syncData();

        Txx->derS.calcDerivative1_x(A11);
        Txy->derS.calcDerivative1_x(A12);
        Tzx->derS.calcDerivative1_x(A13);
        Txy->derS.calcDerivative1_y(A21);
        Tyy->derS.calcDerivative1_y(A22);
        Tyz->derS.calcDerivative1_y(A23);
        Tzx->derS.calcDerivative1_z(A31);
        Tyz->derS.calcDerivative1_z(A32);
        Tzz->derS.calcDerivative1_z(A33);

        A11 = A11 + A21 + A31;
        A22 = A12 + A22 + A32;
        A33 = A13 + A23 + A33;

        xS = V.Vx.fCore.lbound(0) + 2; xE = V.Vx.fCore.ubound(0) - 2;
        yS = V.Vx.fCore.lbound(1) + 2; yE = V.Vx.fCore.ubound(1) - 2;
        zS = V.Vx.fCore.lbound(2) + 2; zE = V.Vx.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vx(iX, iY, iZ) += (A11(iX, iY, iZ) + A11(iX + 1, iY, iZ))*0.5;
                }
            }
        }

        xS = V.Vy.fCore.lbound(0) + 2; xE = V.Vy.fCore.ubound(0) - 2;
        yS = V.Vy.fCore.lbound(1) + 2; yE = V.Vy.fCore.ubound(1) - 2;
        zS = V.Vy.fCore.lbound(2) + 2; zE = V.Vy.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vy(iX, iY, iZ) += (A22(iX, iY, iZ) + A22(iX, iY + 1, iZ))*0.5;
                }
            }
        }

        xS = V.Vz.fCore.lbound(0) + 2; xE = V.Vz.fCore.ubound(0) - 2;
        yS = V.Vz.fCore.lbound(1) + 2; yE = V.Vz.fCore.ubound(1) - 2;
        zS = V.Vz.fCore.lbound(2) + 2; zE = V.Vz.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vz(iX, iY, iZ) += (A33(iX, iY, iZ) + A33(iX, iY, iZ + 1))*0.5;
                }
            }
        }
    }

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
    solveVy(V, nseRHS);
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


/**
 ********************************************************************************************************************************************
 * \brief   Function to advance the solution using Euler method and Implicit Crank-Nicholson method
 *
 *          The non-linear terms are advanced using explicit Euler method, while the duffusion terms are
 *          advanced by semi-implicit Crank-Nicholson method.
 *          This overloaded function advances velocity, temperature and pressure fields for scalar simulations.
 *
 ********************************************************************************************************************************************
 */
void eulerCN_d3::timeAdvance(vfield &V, sfield &P, sfield &T) {
    static plainvf nseRHS(mesh, V);
    static plainsf tmpRHS(mesh, T);

    nseRHS = 0.0;
    tmpRHS = 0.0;

    // First compute the explicit part of the semi-implicit viscous term and divide it by Re
    V.computeDiff(nseRHS);
    nseRHS *= nu;

    // Compute the non-linear term and subtract it from the RHS
    V.computeNLin(V, nseRHS);

    // Add the velocity forcing term
    V.vForcing->addForcing(nseRHS);

    // Add sub-grid stress contribution from LES Model, if enabled
    if (mesh.inputParams.lesModel) {
        // Array limits for loops
        int xS, xE, yS, yE, zS, zE;

        // Temporary variables to store output from spiral LES solver
        double sTxx, sTyy, sTzz, sTxy, sTyz, sTzx;
        double sQx, sQy, sQz;

        // These are three 3x3x3 arrays containing local interpolated velocities
        // These are used to calculate the structure function within the spiral les routine
        blitz::Array<double, 3> u(3, 3, 3), v(3, 3, 3), w(3, 3, 3);

        // The following 3x3 matrix stores the velocity gradient tensor
        blitz::Array<double, 2> dudx(3, 3);

        // The following TinyVector stores the temperature gradient vector
        blitz::TinyVector<double, 3> dtdx;

        // The following TinyVector stores the spiral vortex alignment vector
        blitz::TinyVector<double, 3> e;

        // These 9 arrays store components of the velocity gradient tensor intially
        // Then they are reused to store the derivatives of stress tensor to calculate its divergence
        blitz::Array<double, 3> A11, A12, A13;
        blitz::Array<double, 3> A21, A22, A23;
        blitz::Array<double, 3> A31, A32, A33;

        // These 3 additional arrays are necessary for computing scalar turbulent SGS diffusion
        blitz::Array<double, 3> B1, B2, B3;

        // The 9 arrays of tensor components have the same dimensions and limits as a cell centered variable
        A11.resize(P.F.fSize);      A11.reindexSelf(P.F.flBound);
        A12.resize(P.F.fSize);      A12.reindexSelf(P.F.flBound);
        A13.resize(P.F.fSize);      A13.reindexSelf(P.F.flBound);
        A21.resize(P.F.fSize);      A21.reindexSelf(P.F.flBound);
        A22.resize(P.F.fSize);      A22.reindexSelf(P.F.flBound);
        A23.resize(P.F.fSize);      A23.reindexSelf(P.F.flBound);
        A31.resize(P.F.fSize);      A31.reindexSelf(P.F.flBound);
        A32.resize(P.F.fSize);      A32.reindexSelf(P.F.flBound);
        A33.resize(P.F.fSize);      A33.reindexSelf(P.F.flBound);

        B1.resize(P.F.fSize);       B1.reindexSelf(P.F.flBound);
        B2.resize(P.F.fSize);       B2.reindexSelf(P.F.flBound);
        B3.resize(P.F.fSize);       B3.reindexSelf(P.F.flBound);

        // First interpolate all velocities to cell centers
        // Then U, V and W data are available at all cell centers except ghost points
        P.interTempF = 0.0;
        for (unsigned int i=0; i < P.F.VxIntSlices.size(); i++) {
            P.interTempF(P.F.fCore) += V.Vx.F(P.F.VxIntSlices(i));
        }
        Vxcc->F.F = P.interTempF/P.F.VxIntSlices.size();

        P.interTempF = 0.0;
        for (unsigned int i=0; i < P.F.VyIntSlices.size(); i++) {
            P.interTempF(P.F.fCore) += V.Vy.F(P.F.VyIntSlices(i));
        }
        Vycc->F.F = P.interTempF/P.F.VyIntSlices.size();

        P.interTempF = 0.0;
        for (unsigned int i=0; i < P.F.VzIntSlices.size(); i++) {
            P.interTempF(P.F.fCore) += V.Vz.F(P.F.VzIntSlices(i));
        }
        Vzcc->F.F = P.interTempF/P.F.VzIntSlices.size();

        Vxcc->derS.calcDerivative1_x(A11);
        Vxcc->derS.calcDerivative1_y(A12);
        Vxcc->derS.calcDerivative1_z(A13);
        Vycc->derS.calcDerivative1_x(A21);
        Vycc->derS.calcDerivative1_y(A22);
        Vycc->derS.calcDerivative1_z(A23);
        Vzcc->derS.calcDerivative1_x(A31);
        Vzcc->derS.calcDerivative1_y(A32);
        Vzcc->derS.calcDerivative1_z(A33);

        T.derS.calcDerivative1_x(B1);
        T.derS.calcDerivative1_y(B2);
        T.derS.calcDerivative1_z(B3);

        // Use only cell centers like a collocated grid and compute T tensor
        // Since U, V, and W data is available only in the core,
        // Adjust the loop limits so that the boundary points are excluded
        // to compute derivatives and structure functions correctly.
        xS = P.F.fCore.lbound(0) + 1; xE = P.F.fCore.ubound(0) - 1;
        yS = P.F.fCore.lbound(1) + 1; yE = P.F.fCore.ubound(1) - 1;
        zS = P.F.fCore.lbound(2) + 1; zE = P.F.fCore.ubound(2) - 1;

        for (int iX = xS; iX <= xE; iX++) {
            double dx = mesh.xColloc(iX - 1) - mesh.xColloc(iX);
            for (int iY = yS; iY <= yE; iY++) {
                double dy = mesh.yColloc(iY - 1) - mesh.yColloc(iY);
                for (int iZ = zS; iZ <= zE; iZ++) {
                    double dz = mesh.zColloc(iZ - 1) - mesh.zColloc(iZ);

                    // Cutoff wavelength
                    double del = std::cbrt(dx*dy*dz);

                    u = Vxcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                    v = Vycc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));
                    w = Vzcc->F.F(blitz::Range(iX-1, iX+1), blitz::Range(iY-1, iY+1), blitz::Range(iZ-1, iZ+1));

                    // A 3 x 3 matrix of velocity gradient
                    dudx = A11(iX, iY, iZ), A12(iX, iY, iZ), A13(iX, iY, iZ),
                           A21(iX, iY, iZ), A22(iX, iY, iZ), A23(iX, iY, iZ),
                           A31(iX, iY, iZ), A32(iX, iY, iZ), A33(iX, iY, iZ);

                    double x[3] = {mesh.xStaggr(iX-1), mesh.xStaggr(iX), mesh.xStaggr(iX+1)};
                    double y[3] = {mesh.yStaggr(iY-1), mesh.yStaggr(iY), mesh.yStaggr(iY+1)};
                    double z[3] = {mesh.zStaggr(iZ-1), mesh.zStaggr(iZ), mesh.zStaggr(iZ+1)};

                    sgsLES->sgs_stress(u, v, w, dudx, x, y, z, nu, del,
                                    &sTxx, &sTyy, &sTzz, &sTxy, &sTyz, &sTzx);

                    Txx->F.F(iX, iY, iZ) = sTxx;
                    Tyy->F.F(iX, iY, iZ) = sTyy;
                    Tzz->F.F(iX, iY, iZ) = sTzz;
                    Txy->F.F(iX, iY, iZ) = sTxy;
                    Tyz->F.F(iX, iY, iZ) = sTyz;
                    Tzx->F.F(iX, iY, iZ) = sTzx;

                    // A 3 x 1 vector of temperature gradient
                    dtdx = B1(iX, iY, iZ), B2(iX, iY, iZ), B3(iX, iY, iZ);

                    sgsLES->sgs_flux(dtdx, del, &sQx, &sQy, &sQz);

                    qX->F.F(iX, iY, iZ) = sQx;
                    qY->F.F(iX, iY, iZ) = sQy;
                    qZ->F.F(iX, iY, iZ) = sQz;
                }
            }
        }

        Txx->syncData();
        Tyy->syncData();
        Tzz->syncData();
        Txy->syncData();
        Tyz->syncData();
        Tzx->syncData();

        Txx->derS.calcDerivative1_x(A11);
        Txy->derS.calcDerivative1_x(A12);
        Tzx->derS.calcDerivative1_x(A13);
        Txy->derS.calcDerivative1_y(A21);
        Tyy->derS.calcDerivative1_y(A22);
        Tyz->derS.calcDerivative1_y(A23);
        Tzx->derS.calcDerivative1_z(A31);
        Tyz->derS.calcDerivative1_z(A32);
        Tzz->derS.calcDerivative1_z(A33);

        B1 = A11 + A21 + A31;
        B2 = A12 + A22 + A32;
        B3 = A13 + A23 + A33;

        xS = V.Vx.fCore.lbound(0) + 2; xE = V.Vx.fCore.ubound(0) - 2;
        yS = V.Vx.fCore.lbound(1) + 2; yE = V.Vx.fCore.ubound(1) - 2;
        zS = V.Vx.fCore.lbound(2) + 2; zE = V.Vx.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vx(iX, iY, iZ) += (B1(iX, iY, iZ) + B1(iX + 1, iY, iZ))*0.5;
                }
            }
        }

        xS = V.Vy.fCore.lbound(0) + 2; xE = V.Vy.fCore.ubound(0) - 2;
        yS = V.Vy.fCore.lbound(1) + 2; yE = V.Vy.fCore.ubound(1) - 2;
        zS = V.Vy.fCore.lbound(2) + 2; zE = V.Vy.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vy(iX, iY, iZ) += (B2(iX, iY, iZ) + B2(iX, iY + 1, iZ))*0.5;
                }
            }
        }

        xS = V.Vz.fCore.lbound(0) + 2; xE = V.Vz.fCore.ubound(0) - 2;
        yS = V.Vz.fCore.lbound(1) + 2; yE = V.Vz.fCore.ubound(1) - 2;
        zS = V.Vz.fCore.lbound(2) + 2; zE = V.Vz.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vz(iX, iY, iZ) += (B3(iX, iY, iZ) + B3(iX, iY, iZ + 1))*0.5;
                }
            }
        }

        qX->syncData();
        qY->syncData();
        qZ->syncData();

        qX->derS.calcDerivative1_x(B1);
        qY->derS.calcDerivative1_y(B2);
        qZ->derS.calcDerivative1_z(B3);

        xS = T.F.fCore.lbound(0) + 2; xE = T.F.fCore.ubound(0) - 2;
        yS = T.F.fCore.lbound(1) + 2; yE = T.F.fCore.ubound(1) - 2;
        zS = T.F.fCore.lbound(2) + 2; zE = T.F.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    tmpRHS.F(iX, iY, iZ) = (B1(iX, iY, iX) + B2(iX, iY, iX) + B3(iX, iY, iX));
                }
            }
        }
    }

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
    solveVy(V, nseRHS);
    solveVz(V, nseRHS);

    // Calculate the rhs for the poisson solver (mgRHS) using the divergence of guessed velocity in V
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;

    // Using the calculated mgRHS, evaluate pressure correction (Pp) using multi-grid method
    mgSolver.mgSolve(Pp, mgRHS);

    // Synchronise the pressure correction term across processors
    Pp.syncData();

    // Add the pressure correction term to the pressure field of previous time-step, P
    P += Pp;

    // Finally get the velocity field at end of time-step by subtracting the gradient of pressure correction from V
    Pp.gradient(pressureGradient, V);
    pressureGradient *= dt;
    V -= pressureGradient;

    // Next, for temperature, again compute semi-implicit diffusion term first
    T.computeDiff(tmpRHS);
    tmpRHS *= kappa;

    // Compute the non-linear term and subtract it from the RHS
    T.computeNLin(V, tmpRHS);

    // Add the scalar forcing term
    T.tForcing->addForcing(tmpRHS);

    // Multiply the entire RHS with dt and add the temperature of previous time-step to advance by explicit Euler method
    tmpRHS *= dt;
    tmpRHS += T;

    // Synchronize the RHS term across all processors by updating its sub-domain pads
    tmpRHS.syncData();

    // Using the RHS term computed, compute the guessed temperature of CN method iteratively (and store it in T)
    solveT(T, tmpRHS);

    // Impose boundary conditions on the updated velocity field, V
    V.imposeBCs();

    // Impose boundary conditions on the updated temperature field, T
    T.imposeBCs();
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for x-velocity
 *
 *          The implicit equation for \f$ u_x' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 ********************************************************************************************************************************************
 */
void eulerCN_d3::solveVx(vfield &V, plainvf &nseRHS) {
    int iterCount = 0;
    real locMax = 0.0;
    real gloMax = 0.0;
    static blitz::Array<real, 3> tempVx(V.Vx.F.lbound(), V.Vx.F.shape());

    while (true) {
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(V) shared(nseRHS) shared(tempVx)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    tempVx(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                           hz2hx2 * mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) + V.Vx.F(iX, iY-1, iZ)) +
                                           hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vx(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Colloc(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Staggr(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vx.F = tempVx;

        V.imposeVxBC();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(V) shared(tempVx)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    tempVx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx2) +
                              mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY-1, iZ)) / (hy2) +
                              mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz2));
                }
            }
        }

        tempVx(V.Vx.fBulk) = abs(tempVx(V.Vx.fBulk) - nseRHS.Vx(V.Vx.fBulk));

        locMax = blitz::max(tempVx);
        MPI_Allreduce(&locMax, &gloMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
        if (gloMax < mesh.inputParams.cnTolerance) {
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


/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for y-velocity
 *
 *          The implicit equation for \f$ u_y' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 ********************************************************************************************************************************************
 */
void eulerCN_d3::solveVy(vfield &V, plainvf &nseRHS) {
    int iterCount = 0;
    real locMax = 0.0;
    real gloMax = 0.0;
    static blitz::Array<real, 3> tempVy(V.Vy.F.lbound(), V.Vy.F.shape());

    while (true) {
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(V) shared(nseRHS) shared(tempVy)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    tempVy(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) +
                                           hz2hx2 * mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) + V.Vy.F(iX, iY-1, iZ)) +
                                           hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) + V.Vy.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vy(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Colloc(iY) + hx2hy2 * mesh.ztz2Staggr(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vy.F = tempVy;

        V.imposeVyBC();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(V) shared(tempVy)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    tempVy(iX, iY, iZ) = V.Vy.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) / (hx2) +
                              mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY-1, iZ)) / (hy2) +
                              mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY, iZ-1)) / (hz2));
                }
            }
        }

        tempVy(V.Vy.fBulk) = abs(tempVy(V.Vy.fBulk) - nseRHS.Vy(V.Vy.fBulk));

        locMax = blitz::max(tempVy);
        MPI_Allreduce(&locMax, &gloMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
        if (gloMax < mesh.inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vy not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for z-velocity
 *
 *          The implicit equation for \f$ u_z' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 ********************************************************************************************************************************************
 */
void eulerCN_d3::solveVz(vfield &V, plainvf &nseRHS) {
    int iterCount = 0;
    real locMax = 0.0;
    real gloMax = 0.0;
    static blitz::Array<real, 3> tempVz(V.Vz.F.lbound(), V.Vz.F.shape());

    while (true) {
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(V) shared(nseRHS) shared(tempVz)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    tempVz(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                           hz2hx2 * mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) + V.Vz.F(iX, iY-1, iZ)) +
                                           hx2hy2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vz(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Colloc(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vz.F = tempVz;

        V.imposeVzBC();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(V) shared(tempVz)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    tempVz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx2) +
                              mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY-1, iZ)) / (hy2) +
                              mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz2));
                }
            }
        }

        tempVz(V.Vz.fBulk) = abs(tempVz(V.Vz.fBulk) - nseRHS.Vz(V.Vz.fBulk));

        locMax = blitz::max(tempVz);
        MPI_Allreduce(&locMax, &gloMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
        if (gloMax < mesh.inputParams.cnTolerance) {
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


void eulerCN_d3::solveT(sfield &T, plainsf &tmpRHS) {
    int iterCount = 0;
    real locMax = 0.0;
    real gloMax = 0.0;
    static blitz::Array<real, 3> tempT(T.F.F.lbound(), T.F.F.shape());

    while (true) {
#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(T) shared(tmpRHS) shared(tempT)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
                for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                    tempT(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) + T.F.F(iX-1, iY, iZ)) +
                                          hz2hx2 * mesh.ety2Staggr(iY) * (T.F.F(iX, iY+1, iZ) + T.F.F(iX, iY-1, iZ)) +
                                          hx2hy2 * mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) + T.F.F(iX, iY, iZ-1))) *
                        dt * kappa / ( hx2hy2hz2 * 2.0) + tmpRHS.F(iX, iY, iZ)) /
                 (1.0 + dt * kappa * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Colloc(iZ)))/hx2hy2hz2);
                }
            }
        }

        T.F.F = tempT;

        T.imposeBCs();

#pragma omp parallel for num_threads(mesh.inputParams.nThreads) default(none) shared(T) shared(tempT)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
                for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                    tempT(iX, iY, iZ) = T.F.F(iX, iY, iZ) - 0.5 * dt * kappa * (
                           mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX-1, iY, iZ)) / (hx2) +
                           mesh.ety2Staggr(iY) * (T.F.F(iX, iY+1, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY-1, iZ)) / (hy2) +
                           mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY, iZ-1)) / (hz2));
                }
            }
        }

        tempT(T.F.fBulk) = abs(tempT(T.F.fBulk) - tmpRHS.F(T.F.fBulk));

        locMax = blitz::max(tempT);
        MPI_Allreduce(&locMax, &gloMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
        if (gloMax < mesh.inputParams.cnTolerance) {
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
void eulerCN_d3::setCoefficients() {
    hx2 = pow(mesh.dXi, 2.0);
    hy2 = pow(mesh.dEt, 2.0);
    hz2 = pow(mesh.dZt, 2.0);

    hz2hx2 = pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);
    hx2hy2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0);
    hy2hz2 = pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);

    hx2hy2hz2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
};
