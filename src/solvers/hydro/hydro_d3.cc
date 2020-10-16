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
/*! \file hydro_d3.cc
 *
 *  \brief Definitions of functions for 3D computations with the hydro class.
 *  \sa hydro.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <sys/time.h>
#include <ctime>

#include "hydro.h"
#include "initial.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the hydro_d3 class derived from the base hydro class
 *
 *          The constructor passes its arguments to the base hydro class and then initializes all the scalar and
 *          vector fields necessary for solving the NS equations.
 *          The various coefficients for solving the equations are also set by a call to the \ref setCoefficients function.
 *          Based on the problem type specified by the user in the parameters file, and stored by the \ref parser class as
 *          \ref parser#probType "probType", the appropriate boundary conditions are specified.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
hydro_d3::hydro_d3(const grid &mesh, const parser &solParam, parallel &mpiParam):
            hydro(mesh, solParam, mpiParam),
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
        time = 0.0;

        // INITIALIZE PRESSURE TO 1.0 THROUGHOUT THE DOMAIN
        P = 1.0;

        // INITIALIZE VARIABLES
        initial *initCond;
        switch (inputParams.icType) {
            case 0: initCond = new zeroInitial(mesh);
                break;
            case 1: initCond = new taylorGreen(mesh);
                break;
            case 2: initCond = new channelSine(mesh);
                break;
            case 3: initCond = new uniformRandom(mesh);
                break;
            case 4: initCond = new parabolicRandom(mesh);
                break;
            case 5: initCond = new sineRandom(mesh);
                break;
            default: initCond = new zeroInitial(mesh);
        }
        initCond->initializeField(V);
    }

    checkPeriodic();

    // Initialize velocity boundary conditions
    initVBCs();

    // Initialize velocity forcing field
    initVForcing();

    // Impose boundary conditions on velocity field
    V.imposeBCs();

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
    }
}


void hydro_d3::solvePDE() {
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

    timeStepCount = 0;

    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    tsWriter.writeTSData();

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

        if (timeStepCount % inputParams.ioCnt == 0) {
            tsWriter.writeTSData();
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

        if (std::abs(inputParams.tMax - time) < 0.5*dt) {
            break;
        }
    }
}


void hydro_d3::timeAdvance() {
    nseRHS = 0.0;

    // First compute the explicit part of the semi-implicit viscous term and divide it by Re
    V.computeDiff(nseRHS);
    nseRHS *= inverseRe;

    // Compute the non-linear term and subtract it from the RHS
    V.computeNLin(V, nseRHS);

    // Add the velocity forcing term
    V.vForcing->addForcing(nseRHS);

    // Add sub-grid stress contribution from LES Model, if enabled
    MPI_Barrier(MPI_COMM_WORLD);
    if (inputParams.lesModel and time > 10*inputParams.tStp) {
        // Kinematic viscosity
        double nu = inverseRe;

        // Number of velocity samples for structure function
        int n = 2;

        double v1, v2;

        //if (mesh.rankData.rank == 1) {
        //    std::cout << P.F.fCore.ubound() << std::endl;
        //    std::cout << V.Vx.F.ubound() << std::endl;
        //    std::cout << V.Vy.F.ubound() << std::endl;
        //    std::cout << V.Vz.F.ubound() << std::endl;
        //}
        //MPI_Finalize();
        //exit(0);
        for (int iX = P.F.fCore.lbound(0); iX < P.F.fCore.ubound(0); iX++) {
            double dx = mesh.xColloc(iX - 1) - mesh.xColloc(iX);
            for (int iY = P.F.fCore.lbound(1); iY < P.F.fCore.ubound(1); iY++) {
                double dy = mesh.yColloc(iY - 1) - mesh.yColloc(iY);
                for (int iZ = P.F.fCore.lbound(2); iZ < P.F.fCore.ubound(2); iZ++) {
                    double dz = mesh.zColloc(iZ - 1) - mesh.zColloc(iZ);

                    // Cutoff wavelength
                    double del = std::cbrt(dx*dy*dz);

                    // Vector of alignment of subgrid vortex
                    double e[3] = {0.0, 0.0, 0.0};

                    v1 = (V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))*0.5;
                    v2 = (V.Vx.F(iX, iY+1, iZ+1) + V.Vx.F(iX+1, iY+1, iZ+1))*0.5;
                    double u[2] = {v1, v2};

                    v1 = (V.Vy.F(iX, iY-1, iZ) + V.Vy.F(iX, iY, iZ))*0.5;
                    v2 = (V.Vy.F(iX+1, iY, iZ+1) + V.Vy.F(iX+1, iY+1, iZ+1))*0.5;
                    double v[2] = {v1, v2};

                    v1 = (V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))*0.5;
                    v2 = (V.Vz.F(iX+1, iY+1, iZ) + V.Vz.F(iX+1, iY+1, iZ+1))*0.5;
                    double w[2] = {v1, v2};

                    double ux = (V.Vx.F(iX, iY, iZ) - V.Vx.F(iX-1, iY, iZ))/(mesh.xColloc(iX) - mesh.xColloc(iX-1));
                    double vy = (V.Vy.F(iX, iY, iZ) - V.Vy.F(iX, iY-1, iZ))/(mesh.yColloc(iY) - mesh.yColloc(iY-1));
                    double wz = (V.Vz.F(iX, iY, iZ) - V.Vz.F(iX, iY, iZ-1))/(mesh.zColloc(iZ) - mesh.zColloc(iZ-1));

                    v1 = (V.Vx.F(iX-1, iY-1, iZ) + V.Vx.F(iX, iY-1, iZ))*0.5;
                    v2 = (V.Vx.F(iX-1, iY+1, iZ) + V.Vx.F(iX, iY+1, iZ))*0.5;
                    double uy = (v2 - v1)/(mesh.yStaggr(iY+1) - mesh.yStaggr(iY-1));

                    v1 = (V.Vx.F(iX-1, iY, iZ-1) + V.Vx.F(iX, iY, iZ-1))*0.5;
                    v2 = (V.Vx.F(iX-1, iY, iZ+1) + V.Vx.F(iX, iY, iZ+1))*0.5;
                    double uz = (v2 - v1)/(mesh.zStaggr(iZ+1) - mesh.zStaggr(iZ-1));

                    v1 = (V.Vy.F(iX-1, iY-1, iZ) + V.Vy.F(iX-1, iY, iZ))*0.5;
                    v2 = (V.Vy.F(iX+1, iY-1, iZ) + V.Vy.F(iX+1, iY, iZ))*0.5;
                    double vx = (v2 - v1)/(mesh.xStaggr(iX+1) - mesh.xStaggr(iX-1));

                    v1 = (V.Vy.F(iX, iY-1, iZ-1) + V.Vy.F(iX, iY, iZ-1))*0.5;
                    v2 = (V.Vy.F(iX, iY-1, iZ+1) + V.Vy.F(iX, iY, iZ+1))*0.5;
                    double vz = (v2 - v1)/(mesh.zStaggr(iZ+1) - mesh.zStaggr(iZ-1));

                    v1 = (V.Vz.F(iX-1, iY, iZ-1) + V.Vz.F(iX-1, iY, iZ))*0.5;
                    v2 = (V.Vz.F(iX+1, iY, iZ-1) + V.Vz.F(iX+1, iY, iZ))*0.5;
                    double wx = (v2 - v1)/(mesh.xStaggr(iX+1) - mesh.xStaggr(iX-1));

                    v1 = (V.Vz.F(iX, iY-1, iZ-1) + V.Vz.F(iX, iY-1, iZ))*0.5;
                    v2 = (V.Vz.F(iX, iY+1, iZ-1) + V.Vz.F(iX, iY+1, iZ))*0.5;
                    double wy = (v2 - v1)/(mesh.yStaggr(iY+1) - mesh.yStaggr(iY-1));

                    // A 3 x 3 matrix of velocity gradient
                    double dudx[3][3] = {{ux, uy, uz}, {vx, vy, vz}, {wx, wy, wz}};

                    double x[2] = {mesh.xStaggr(iX), mesh.xStaggr(iX+1)};
                    double y[2] = {mesh.yStaggr(iY), mesh.yStaggr(iY+1)};
                    double z[2] = {mesh.zStaggr(iZ), mesh.zStaggr(iZ+1)};

                    double sTxx, sTyy, sTzz, sTxy, sTyz, sTzx;

                    sgsLES->sgs_stress(u, v, w, x, y, z, n, dudx, e, nu, del,
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

        blitz::Array<double, 3> dTxxDx;
        dTxxDx.resize(Txx->F.fSize);
        dTxxDx.reindexSelf(Txx->F.flBound);
        Txx->derS.calcDerivative1_x(dTxxDx);

        blitz::Array<double, 3> dTxyDx;
        dTxyDx.resize(Txy->F.fSize);
        dTxyDx.reindexSelf(Txy->F.flBound);
        Txy->derS.calcDerivative1_x(dTxyDx);

        blitz::Array<double, 3> dTzxDx;
        dTzxDx.resize(Tzx->F.fSize);
        dTzxDx.reindexSelf(Tzx->F.flBound);
        Tzx->derS.calcDerivative1_x(dTzxDx);

        blitz::Array<double, 3> dTxyDy;
        dTxyDy.resize(Txy->F.fSize);
        dTxyDy.reindexSelf(Txy->F.flBound);
        Txy->derS.calcDerivative1_y(dTxyDy);

        blitz::Array<double, 3> dTyyDy;
        dTyyDy.resize(Tyy->F.fSize);
        dTyyDy.reindexSelf(Tyy->F.flBound);
        Tyy->derS.calcDerivative1_y(dTyyDy);

        blitz::Array<double, 3> dTyzDy;
        dTyzDy.resize(Tyz->F.fSize);
        dTyzDy.reindexSelf(Tyz->F.flBound);
        Tyz->derS.calcDerivative1_y(dTyzDy);

        blitz::Array<double, 3> dTzxDz;
        dTzxDz.resize(Tzx->F.fSize);
        dTzxDz.reindexSelf(Tzx->F.flBound);
        Tzx->derS.calcDerivative1_z(dTzxDz);

        blitz::Array<double, 3> dTyzDz;
        dTyzDz.resize(Tyz->F.fSize);
        dTyzDz.reindexSelf(Tyz->F.flBound);
        Tyz->derS.calcDerivative1_z(dTyzDz);

        blitz::Array<double, 3> dTzzDz;
        dTzzDz.resize(Tzz->F.fSize);
        dTzzDz.reindexSelf(Tzz->F.flBound);
        Tzz->derS.calcDerivative1_z(dTzzDz);

        dTxxDx = dTxxDx + dTxyDy + dTzxDz;
        dTyyDy = dTxyDx + dTyyDy + dTyzDz;
        dTzzDz = dTzxDx + dTyzDy + dTzzDz;

        int xS, xE, yS, yE, zS, zE;

        xS = V.Vx.fCore.lbound(0) + 2; xE = V.Vx.fCore.ubound(0) - 2;
        yS = V.Vx.fCore.lbound(1) + 2; yE = V.Vx.fCore.ubound(1) - 2;
        zS = V.Vx.fCore.lbound(2) + 2; zE = V.Vx.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vx(iX, iY, iZ) += (dTxxDx(iX, iY, iZ) + dTxxDx(iX + 1, iY, iZ))*0.5;
                }
            }
        }

        xS = V.Vy.fCore.lbound(0) + 2; xE = V.Vy.fCore.ubound(0) - 2;
        yS = V.Vy.fCore.lbound(1) + 2; yE = V.Vy.fCore.ubound(1) - 2;
        zS = V.Vy.fCore.lbound(2) + 2; zE = V.Vy.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vy(iX, iY, iZ) += (dTyyDy(iX, iY, iZ) + dTyyDy(iX, iY + 1, iZ))*0.5;
                }
            }
        }

        xS = V.Vz.fCore.lbound(0) + 2; xE = V.Vz.fCore.ubound(0) - 2;
        yS = V.Vz.fCore.lbound(1) + 2; yE = V.Vz.fCore.ubound(1) - 2;
        zS = V.Vz.fCore.lbound(2) + 2; zE = V.Vz.fCore.ubound(2) - 2;

        for (int iX = xS; iX <= xE; iX++) {
            for (int iY = yS; iY <= yE; iY++) {
                for (int iZ = zS; iZ <= zE; iZ++) {
                    nseRHS.Vz(iX, iY, iZ) += (dTzzDz(iX, iY, iZ) + dTzzDz(iX, iY, iZ + 1))*0.5;
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
    solveVx();
    solveVy();
    solveVz();

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

void hydro_d3::solveVx() {
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
                                     dt / ( hx2hy2hz2 * 2.0 * inputParams.Re) + nseRHS.Vx(iX, iY, iZ)) /
                                 (1.0 + dt * ((hy2hz2 * mesh.xix2Colloc(iX) +
                                               hz2hx2 * mesh.ety2Staggr(iY) +
                                               hx2hy2 * mesh.ztz2Staggr(iZ)))/(inputParams.Re * hx2hy2hz2));
                }
            }
        }

        V.Vx.F = guessedVelocity.Vx;

        V.imposeVxBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - 0.5 * dt * inverseRe * (
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

void hydro_d3::solveVy() {
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
                                     dt / ( hx2hy2hz2 * 2.0 * inputParams.Re) + nseRHS.Vy(iX, iY, iZ)) /
                                 (1.0 + dt * ((hy2hz2 * mesh.xix2Staggr(iX) +
                                               hz2hx2 * mesh.ety2Colloc(iY) +
                                               hx2hy2 * mesh.ztz2Staggr(iZ)))/(inputParams.Re * hx2hy2hz2));
                }
            }
        }

        V.Vy.F = guessedVelocity.Vy;

        V.imposeVyBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vy(iX, iY, iZ) = V.Vy.F(iX, iY, iZ) - 0.5 * dt * inverseRe * (
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

void hydro_d3::solveVz() {
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
                                     dt / ( hx2hy2hz2 * 2.0 * inputParams.Re) + nseRHS.Vz(iX, iY, iZ)) /
                                 (1.0 + dt * ((hy2hz2 * mesh.xix2Staggr(iX) +
                                               hz2hx2 * mesh.ety2Staggr(iY) +
                                               hx2hy2 * mesh.ztz2Colloc(iZ)))/(inputParams.Re * hx2hy2hz2));
                }
            }
        }

        V.Vz.F = guessedVelocity.Vz;

        V.imposeVzBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - 0.5 * dt * inverseRe * (
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


real hydro_d3::testPeriodic() {
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

hydro_d3::~hydro_d3() { }
