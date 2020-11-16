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
/*! \file scalar_d3.cc
 *
 *  \brief Definitions of functions for 3D computations with the scalar class.
 *  \sa scalar.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <sys/time.h>
#include <ctime>

#include "scalar.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the scalar_d3 class derived from the base scalar class
 *
 *          The constructor passes its arguments to the base scalar class and then initializes all the scalar and
 *          vector fields necessary for solving the NS equations.
 *          The various coefficients for solving the equations are also set by a call to the \ref setCoefficients function.
 *          Based on the problem type specified by the user in the parameters file, and stored by the \ref parser class as
 *          \ref parser#probType "probType", the appropriate boundary conditions are specified.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
scalar_d3::scalar_d3(const grid &mesh, const parser &solParam, parallel &mpiParam):
            scalar(mesh, solParam, mpiParam),
            mgSolver(mesh, inputParams)
{
    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // INITIALIZE VARIABLES
    if (inputParams.restartFlag) {
        // Fields to be read from HDF5 file are passed to reader class as a vector
        std::vector<field> readFields;

        // Populate the vector with required fields
        readFields.push_back(V.Vx);
        readFields.push_back(V.Vy);
        readFields.push_back(V.Vz);
        readFields.push_back(P.F);
        readFields.push_back(T.F);

        // Initialize reader object
        reader dataReader(mesh, readFields);

        time = dataReader.readData();

        // Abort if this time is greater than the final time specified by the user
        if (time >= inputParams.tMax) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Restart file is starting from a point beyond the final time specified. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    } else {
        P = 1.0;
        T = 0.0;
        time = 0.0;
    }

    checkPeriodic();

    // Initialize velocity and temperature boundary conditions
    initVBCs();
    initTBCs();

    // Initialize velocity and temperature forcing fields
    initVForcing();
    initTForcing();

    // Impose boundary conditions on velocity and temperature fields
    V.imposeBCs();
    T.imposeBCs();

    // If LES switch is enabled, initialize LES model
    if (inputParams.lesModel) {
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


void scalar_d3::solvePDE() {
#ifdef TIME_RUN
    visc_time = 0.0;
    nlin_time = 0.0;
    intr_time = 0.0;
    impl_time = 0.0;
    prhs_time = 0.0;
    pois_time = 0.0;

#else
    real fwTime, prTime, rsTime;

    // Set dt equal to input time step
    dt = inputParams.tStp;

    // Fields to be written into HDF5 file are passed to writer class as a vector
    std::vector<field> writeFields;

    // Populate the vector with required fields
    writeFields.push_back(V.Vx);
    writeFields.push_back(V.Vy);
    writeFields.push_back(V.Vz);
    writeFields.push_back(P.F);
    writeFields.push_back(T.F);

    // Initialize writer object
    writer dataWriter(mesh, writeFields);

    // Initialize probes
    if (inputParams.readProbes) {
        dataProbe = new probes(mesh, writeFields);
    }

    // Output file containing the time series of various variables
    tseries tsWriter(mesh, V, P, time, dt);

    // FILE WRITING TIME
    fwTime = time;

    // FIELD PROBING TIME
    prTime = time;

    // RESTART FILE WRITING TIME
    rsTime = time;
#endif

    timeStepCount = 0;

#ifndef TIME_RUN
    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    tsWriter.writeTSData(T, nu, kappa);

    // WRITE DATA AT t = 0 OR INCREMENT INTERVAL IF RESTARTING
    if (inputParams.restartFlag) {
        int tCount, fCount;

        tCount = int(time/inputParams.tStp);

        fCount = int(inputParams.fwInt/inputParams.tStp);
        fwTime = roundNum(tCount, fCount)*inputParams.tStp;

        fCount = int(inputParams.prInt/inputParams.tStp);
        prTime = roundNum(tCount, fCount)*inputParams.tStp;

        fCount = int(inputParams.rsInt/inputParams.tStp);
        rsTime = roundNum(tCount, fCount)*inputParams.tStp;
    }

    switch (inputParams.solnFormat) {
        case 1: dataWriter.writeSolution(time);
            break;
        case 2: dataWriter.writeTarang(time);
            break;
        default: dataWriter.writeSolution(time);
    }
    fwTime += inputParams.fwInt;

    if (inputParams.readProbes) {
        dataProbe->probeData(time);
        prTime += inputParams.prInt;
    }

    rsTime += inputParams.rsInt;
#endif

    // TIME-INTEGRATION LOOP
    while (true) {
        // MAIN FUNCTION CALLED IN EACH LOOP TO UPDATE THE FIELDS AT EACH TIME-STEP
        timeAdvance();

        if (inputParams.useCFL) {
            V.computeTStp(dt);
            if (dt > inputParams.tStp) {
                dt = inputParams.tStp;
            }
        }

        timeStepCount += 1;
        time += dt;

#ifndef TIME_RUN
        if (timeStepCount % inputParams.ioCnt == 0) {
            tsWriter.writeTSData(T, nu, kappa);
        }

        if (inputParams.readProbes and std::abs(prTime - time) < 0.5*dt) {
            dataProbe->probeData(time);
            prTime += inputParams.prInt;
        }

        if (std::abs(fwTime - time) < 0.5*dt) {
            switch (inputParams.solnFormat) {
                case 1: dataWriter.writeSolution(time);
                    break;
                case 2: dataWriter.writeTarang(time);
                    break;
                default: dataWriter.writeSolution(time);
            }
            fwTime += inputParams.fwInt;
        }

        if (std::abs(rsTime - time) < 0.5*dt) {
            dataWriter.writeRestart(time);
            rsTime += inputParams.rsInt;
        }
#endif

        if (std::abs(inputParams.tMax - time) < 0.5*dt) {
            break;
        }
    }

    // WRITE THE OUTPUT OF THE TIMING RUN
#ifdef TIME_RUN
    if (mesh.rankData.rank == 0) {
        std::cout << std::left << std::setw(50) << "Time taken in computing viscous terms: "            << std::fixed << std::setprecision(6) << visc_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing non-linear terms: "         << std::fixed << std::setprecision(6) << nlin_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing intermediate steps: "       << std::fixed << std::setprecision(6) << intr_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing velocities implicitly: "    << std::fixed << std::setprecision(6) << impl_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing RHS for poisson solver: "   << std::fixed << std::setprecision(6) << prhs_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken by poisson solver: "                     << std::fixed << std::setprecision(6) << pois_time << std::endl;
    }
#endif
}


void scalar_d3::timeAdvance() {
#ifdef TIME_RUN
    struct timeval begin, end;
#endif

    nseRHS = 0.0;
    tmpRHS = 0.0;

    // First compute the explicit part of the semi-implicit viscous term and divide it by Re
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    V.computeDiff(nseRHS);
    nseRHS *= nu;
    gettimeofday(&end, NULL);
    visc_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    V.computeDiff(nseRHS);
    nseRHS *= nu;
#endif

    // Compute the non-linear term and subtract it from the RHS
#ifndef TIME_RUN
    V.computeNLin(V, nseRHS);
#else
    gettimeofday(&begin, NULL);
    V.computeNLin(V, nseRHS);
    gettimeofday(&end, NULL);
    nlin_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

    gettimeofday(&begin, NULL);
#endif

    // Add the velocity forcing term
    V.vForcing->addForcing(nseRHS);

    // Add sub-grid stress contribution from LES Model, if enabled
    if (inputParams.lesModel and time > 10*inputParams.tStp) {
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

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    intr_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif

    // Synchronize the RHS term across all processors by updating its sub-domain pads
    nseRHS.syncData();

    // Using the RHS term computed, compute the guessed velocity of CN method iteratively (and store it in V)
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
#endif
    solveVx();
    solveVy();
    solveVz();

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    impl_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif

    // Calculate the rhs for the poisson solver (mgRHS) using the divergence of guessed velocity in V
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;
    gettimeofday(&end, NULL);
    prhs_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;
#endif

    // Using the calculated mgRHS, evaluate pressure correction (Pp) using multi-grid method
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    mgSolver.mgSolve(Pp, mgRHS);
    gettimeofday(&end, NULL);
    pois_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    mgSolver.mgSolve(Pp, mgRHS);
#endif

    // Synchronise the pressure correction term across processors
    Pp.syncData();

    // Add the pressure correction term to the pressure field of previous time-step, P
    P += Pp;

    // Finally get the velocity field at end of time-step by subtracting the gradient of pressure correction from V
    Pp.gradient(pressureGradient, V);
    pressureGradient *= dt;
    V -= pressureGradient;

    // Next compute the explicit part of the semi-implicit diffusion term and multiply it by kappa
    T.computeDiff(tmpRHS);
    tmpRHS *= kappa;

    // COMPUTE THE CONVECTIVE DERIVATIVE AND SUBTRACT IT FROM THE CALCULATED DIFFUSION TERMS OF RHS IN tmpRHS
    T.computeNLin(V, tmpRHS);

    // Add the scalar forcing term
    T.tForcing->addForcing(tmpRHS);

    // Multiply the entire RHS with dt and add the temperature of previous time-step to advance by explicit Euler method
    tmpRHS *= dt;
    tmpRHS += T;

    // Synchronize the RHS term across all processors by updating its sub-domain pads
    tmpRHS.syncData();

    // Using the RHS term computed, compute the guessed temperature of CN method iteratively (and store it in T)
    solveT();

    // Impose boundary conditions on the updated velocity field, V
    V.imposeBCs();

    // Impose boundary conditions on the updated temperature field, T
    T.imposeBCs();
}


void scalar_d3::solveVx() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vx(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                                       hz2hx2 * mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) + V.Vx.F(iX, iY-1, iZ)) +
                                                       hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vx(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Colloc(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Staggr(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vx.F = guessedVelocity.Vx;

        V.imposeVxBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx * hx) +
                              mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY-1, iZ)) / (hy * hy) +
                              mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        velocityLaplacian.Vx(V.Vx.fBulk) = abs(velocityLaplacian.Vx(V.Vx.fBulk) - nseRHS.Vx(V.Vx.fBulk));

        maxError = velocityLaplacian.vxMax();

        if (maxError < inputParams.cnTolerance) {
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


void scalar_d3::solveVy() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vy(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) +
                                                       hz2hx2 * mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) + V.Vy.F(iX, iY-1, iZ)) +
                                                       hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) + V.Vy.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vy(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Colloc(iY) + hx2hy2 * mesh.ztz2Staggr(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vy.F = guessedVelocity.Vy;

        V.imposeVyBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vy(iX, iY, iZ) = V.Vy.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) / (hx * hx) +
                              mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY-1, iZ)) / (hy * hy) +
                              mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        velocityLaplacian.Vy(V.Vy.fBulk) = abs(velocityLaplacian.Vy(V.Vy.fBulk) - nseRHS.Vy(V.Vy.fBulk));

        maxError = velocityLaplacian.vyMax();

        if (maxError < inputParams.cnTolerance) {
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


void scalar_d3::solveVz() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vz(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                                       hz2hx2 * mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) + V.Vz.F(iX, iY-1, iZ)) +
                                                       hx2hy2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vz(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Colloc(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vz.F = guessedVelocity.Vz;

        V.imposeVzBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx * hx) +
                              mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY-1, iZ)) / (hy * hy) +
                              mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        velocityLaplacian.Vz(V.Vz.fBulk) = abs(velocityLaplacian.Vz(V.Vz.fBulk) - nseRHS.Vz(V.Vz.fBulk));

        maxError = velocityLaplacian.vzMax();

        if (maxError < inputParams.cnTolerance) {
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


void scalar_d3::solveT() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
                for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                    guessedScalar.F(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) + T.F.F(iX-1, iY, iZ)) +
                                                    hz2hx2 * mesh.ety2Staggr(iY) * (T.F.F(iX, iY+1, iZ) + T.F.F(iX, iY-1, iZ)) +
                                                    hx2hy2 * mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) + T.F.F(iX, iY, iZ-1))) *
                             dt * kappa / ( hx2hy2hz2 * 2.0) + tmpRHS.F(iX, iY, iZ)) /
                      (1.0 + dt * kappa * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Colloc(iZ)))/hx2hy2hz2);
                }
            }
        }

        T = guessedScalar;

        T.imposeBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
                for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                    scalarLaplacian.F(iX, iY, iZ) = T.F.F(iX, iY, iZ) - 0.5 * dt * kappa * (
                           mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX-1, iY, iZ)) / (hx * hx) +
                           mesh.ety2Staggr(iY) * (T.F.F(iX, iY+1, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY-1, iZ)) / (hy * hy) +
                           mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        scalarLaplacian.F(T.F.fBulk) = abs(scalarLaplacian.F(T.F.fBulk) - tmpRHS.F(T.F.fBulk));

        maxError = scalarLaplacian.fxMax();

        if (maxError < inputParams.cnTolerance) {
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


real scalar_d3::testPeriodic() {
    real xCoord = 0.0;
    real yCoord = 0.0;
    real zCoord = 0.0;

    nseRHS = 0.0;
    V = 0.0;

    for (int i=V.Vx.F.lbound(0); i <= V.Vx.F.ubound(0); i++) {
        for (int j=V.Vx.F.lbound(1); j <= V.Vx.F.ubound(1); j++) {
            for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
                V.Vx.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                  cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                  cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                nseRHS.Vx(i, j, k) = V.Vx.F(i, j, k);
            }
        }
    }

    for (int i=V.Vy.F.lbound(0); i <= V.Vy.F.ubound(0); i++) {
        for (int j=V.Vy.F.lbound(1); j <= V.Vy.F.ubound(1); j++) {
            for (int k=V.Vy.F.lbound(2); k <= V.Vy.F.ubound(2); k++) {
                V.Vy.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                   sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                   cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                nseRHS.Vy(i, j, k) = V.Vy.F(i, j, k);
            }
        }
    }

    V.Vz.F = 0.0;
    nseRHS.Vz = 0.0;

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    // X-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iY = V.Vx.fCore.lbound(1); iY <= V.Vx.fCore.ubound(1); iY += iX) {
            for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += iX) {
                xCoord = mesh.xColloc(V.Vx.fCore.lbound(0)) - (mesh.xColloc(V.Vx.fCore.lbound(0) + iX) - mesh.xColloc(V.Vx.fCore.lbound(0)));
                nseRHS.Vx(V.Vx.fCore.lbound(0) - iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xColloc(V.Vx.fCore.ubound(0)) + (mesh.xColloc(V.Vx.fCore.ubound(0)) - mesh.xColloc(V.Vx.fCore.ubound(0) - iX));
                nseRHS.Vx(V.Vx.fCore.ubound(0) + iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // X-VELOCITY IN FRONT AND BACK PADS
    for (int iY = 1; iY <= mesh.padWidths(1); iY++) {
        for (int iX = V.Vx.fCore.lbound(0); iX <= V.Vx.fCore.ubound(0); iX += iY) {
            for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += iY) {
                yCoord = mesh.yStaggr(V.Vx.fCore.lbound(1)) - (mesh.yStaggr(V.Vx.fCore.lbound(1) + iY) - mesh.yStaggr(V.Vx.fCore.lbound(1)));
                nseRHS.Vx(iX, V.Vx.fCore.lbound(1) - iY, iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yStaggr(V.Vx.fCore.ubound(1)) + (mesh.yStaggr(V.Vx.fCore.ubound(1)) - mesh.yStaggr(V.Vx.fCore.ubound(1) - iY));
                nseRHS.Vx(iX, V.Vx.fCore.ubound(1) + iY, iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // X-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vx.fCore.lbound(0); iX <= V.Vx.fCore.ubound(0); iX += iZ) {
            for (int iY = V.Vx.fCore.lbound(1); iY <= V.Vx.fCore.ubound(1); iY += iZ) {
                zCoord = mesh.zStaggr(V.Vx.fCore.lbound(2)) - (mesh.zStaggr(V.Vx.fCore.lbound(2) + iZ) - mesh.zStaggr(V.Vx.fCore.lbound(2)));
                nseRHS.Vx(iX, iY, V.Vx.fCore.lbound(2) - iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(V.Vx.fCore.ubound(2)) + (mesh.zStaggr(V.Vx.fCore.ubound(2)) - mesh.zStaggr(V.Vx.fCore.ubound(2) - iZ));
                nseRHS.Vx(iX, iY, V.Vx.fCore.ubound(2) + iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iY = V.Vy.fCore.lbound(1); iY <= V.Vy.fCore.ubound(1); iY += iX) {
            for (int iZ = V.Vy.fCore.lbound(2); iZ <= V.Vy.fCore.ubound(2); iZ += iX) {
                xCoord = mesh.xStaggr(V.Vy.fCore.lbound(0)) - (mesh.xStaggr(V.Vy.fCore.lbound(0) + iX) - mesh.xStaggr(V.Vy.fCore.lbound(0)));
                nseRHS.Vy(V.Vy.fCore.lbound(0) - iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xStaggr(V.Vy.fCore.ubound(0)) + (mesh.xStaggr(V.Vy.fCore.ubound(0)) - mesh.xStaggr(V.Vy.fCore.ubound(0) - iX));
                nseRHS.Vy(V.Vy.fCore.ubound(0) + iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY IN FRONT AND BACK PADS
    for (int iY = 1; iY <= mesh.padWidths(1); iY++) {
        for (int iX = V.Vy.fCore.lbound(0); iX <= V.Vy.fCore.ubound(0); iX += iY) {
            for (int iZ = V.Vy.fCore.lbound(2); iZ <= V.Vy.fCore.ubound(2); iZ += iY) {
                yCoord = mesh.yColloc(V.Vy.fCore.lbound(1)) - (mesh.yColloc(V.Vy.fCore.lbound(1) + iY) - mesh.yColloc(V.Vy.fCore.lbound(1)));
                nseRHS.Vy(iX, V.Vy.fCore.lbound(1) - iY, iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*yCoord/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yColloc(V.Vy.fCore.ubound(1)) + (mesh.yColloc(V.Vy.fCore.ubound(1)) - mesh.yColloc(V.Vy.fCore.ubound(1) - iY));
                nseRHS.Vy(iX, V.Vy.fCore.ubound(1) + iY, iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*yCoord/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vy.fCore.lbound(0); iX <= V.Vy.fCore.ubound(0); iX += iZ) {
            for (int iY = V.Vy.fCore.lbound(1); iY <= V.Vy.fCore.ubound(1); iY += iZ) {
                zCoord = mesh.zStaggr(V.Vy.fCore.lbound(2)) - (mesh.zStaggr(V.Vy.fCore.lbound(2) + iZ) - mesh.zStaggr(V.Vy.fCore.lbound(2)));
                nseRHS.Vy(iX, iY, V.Vy.fCore.lbound(2) - iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(V.Vy.fCore.ubound(2)) + (mesh.zStaggr(V.Vy.fCore.ubound(2)) - mesh.zStaggr(V.Vy.fCore.ubound(2) - iZ));
                nseRHS.Vy(iX, iY, V.Vy.fCore.ubound(2) + iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    V.imposeBCs();

    V -= nseRHS;

    return std::max(blitz::max(fabs(V.Vx.F)), blitz::max(fabs(V.Vy.F)));
}


scalar_d3::~scalar_d3() { }
