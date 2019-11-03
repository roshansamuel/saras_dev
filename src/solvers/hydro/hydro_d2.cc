#include "hydro.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the hydro_d2 class derived from the base hydro class
 *
 *          The constructor passes its arguments to the base hydro class and then initializes all the scalar and
 *          vector fields necessary for solving the NS equations.
 *          The various coefficients for solving the equations are also set by a call to the \ref setCoefficients function.
 *          Based on the problem type specified by the user in the parameters file, and stored by the \ref parser class as
 *          \ref parser#probType "probType", the appropriate initial and boundary conditions are specified.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
hydro_d2::hydro_d2(const grid &mesh, const parser &solParam, parallel &mpiParam):
            hydro(mesh, solParam, mpiParam),
            mgSolver(mesh, inputParams)
{
    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // INITIALIZE VARIABLES
    if (inputParams.restartFlag) {
        // Scalar fields to be read from HDF5 file is passed to reader class as a vector
        std::vector<field> readFields;

        // Populate the vector with required scalar fields
        readFields.push_back(V.Vx);
        readFields.push_back(V.Vz);
        readFields.push_back(P.F);

        // Initialize reader object
        reader dataReader(mesh, readFields);

        time = dataReader.readData();

    } else {
        time = 0.0;

        // INITIALIZE VARIABLES
        if (inputParams.probType == 1) {
            // INITIALIZE PRESSURE TO 1.0 THROUGHOUT THE DOMAIN
            P = 1.0;

            // FOR LDC, SET THE X-VELOCITY OF STAGGERED POINTS ON THE LID TO 1.0
            V.Vx.F(blitz::Range::all(), blitz::Range::all(), mesh.staggrCoreDomain.ubound(2)) = 1.0;

        } else if (inputParams.probType == 2) {
            // INITIALIZE PRESSURE TO 1.0 THROUGHOUT THE DOMAIN
            P = 1.0;

            // X-VELOCITY
            for (int i=V.Vx.F.lbound(0); i <= V.Vx.F.ubound(0); i++) {
                for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
                    V.Vx.F(i, 0, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                        cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                }
            }

            // Z-VELOCITY
            for (int i=V.Vz.F.lbound(0); i <= V.Vz.F.ubound(0); i++) {
                for (int k=V.Vz.F.lbound(2); k <= V.Vz.F.ubound(2); k++) {
                    V.Vz.F(i, 0, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                         sin(2.0*M_PI*mesh.zColloc(k)/mesh.zLen);
                }
            }
        }
    }

    // Disable periodic data transfer by setting neighbouring ranks of boundary sub-domains to NULL
    // Left and right walls
    if (not inputParams.xPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along X Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.xRank == 0)             mpiData.nearRanks(0) = MPI_PROC_NULL;
        if (mpiData.xRank == mpiData.npX-1) mpiData.nearRanks(1) = MPI_PROC_NULL;
    }

    // Front and rear walls are by default non-periodic for 2D simulations
    if (mpiData.yRank == 0)             mpiData.nearRanks(2) = MPI_PROC_NULL;
    if (mpiData.yRank == mpiData.npY-1) mpiData.nearRanks(3) = MPI_PROC_NULL;

    // Inform user about BC along top and bottom walls
    if (not inputParams.zPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Z Direction" << std::endl;
            std::cout << std::endl;
        }
    }

    imposeUBCs();
    imposeWBCs();
}

void hydro_d2::solvePDE() {
    int xLow, xTop;
    int zLow, zTop;

    double dVol;
    double maxDivergence;
    double fwTime, prTime;
    double totalEnergy, localEnergy;

    // Scalar field to compute divergence
    sfield divV(mesh, "DIV_V", false);

    // Scalar fields to be written into HDF5 file is passed to writer class as a vector
    std::vector<field> writeFields;

    // Populate the vector with required scalar fields
    writeFields.push_back(V.Vx);
    writeFields.push_back(V.Vz);
    writeFields.push_back(P.F);

    // Initialize writer object
    writer dataWriter(mesh, writeFields);

    // Initialize probes
    if (inputParams.readProbes) {
        dataProbe = new probes(mesh, writeFields);
    }

    // Output file containing the time series of various variables
    std::ofstream ofFile("output/TimeSeries.dat");

    // FILE WRITING TIME
    fwTime = time;

    // FIELD PROBING TIME
    prTime = time;

    timeStepCount = 0;

    // PARAMETERS FOR COMPUTING TOTAL ENERGY IN THE DOMAIN
    // UPPER AND LOWER LIMITS WHEN COMPUTING ENERGY IN STAGGERED GRID
    xLow = P.F.fCore.lbound(0);        xTop = P.F.fCore.ubound(0);
    zLow = P.F.fCore.lbound(2);        zTop = P.F.fCore.ubound(2);

    // STAGGERED GRIDS HAVE SHARED POINT ACROSS MPI-SUBDOMAINS - ACCORDINGLY DECREASE LIMITS
    if (mesh.rankData.xRank > 0) {
        xLow += 1;
    }

    // INFINITESIMAL VOLUME FOR INTEGRATING ENERGY OVER DOMAIN
    dVol = hx*hz;

    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    V.divergence(divV);
    maxDivergence = divV.F.fieldMax();

    int iY = 0;
    localEnergy = 0.0;
    totalEnergy = 0.0;
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iZ = zLow; iZ <= zTop; iZ++) {
            localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                            pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*dVol;
        }
    }

    MPI_Allreduce(&localEnergy, &totalEnergy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);

    if (mesh.rankData.rank == 0) {
        std::cout << std::fixed << std::setw(6)  << "Time" << "\t" <<
                                   std::setw(16) << "Total energy" << "\t" << "Divergence" << std::endl;

        std::cout << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                   std::setw(16) << std::setprecision(8) << sqrt(totalEnergy) << "\t" << maxDivergence << std::endl;

        ofFile << "#VARIABLES = Time, Energy, Divergence, dt\n";
        ofFile << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                std::setw(16) << std::setprecision(8) << sqrt(totalEnergy) << "\t" << maxDivergence << "\t" << inputParams.tStp << std::endl;
    }

    // WRITE DATA AT t = 0
    if (not inputParams.restartFlag) {
        dataWriter.writeData(time);

        if (inputParams.readProbes) {
            dataProbe->probeData(time);
        }
    }
    fwTime += inputParams.fwInt;
    prTime += inputParams.prInt;

    dt = inputParams.tStp;
    // TIME-INTEGRATION LOOP
    while (true) {
        // MAIN FUNCTION CALLED IN EACH LOOP TO UPDATE THE FIELDS AT EACH TIME-STEP
        computeTimeStep();
        V.computeTStp(dt);
        if (dt > inputParams.tStp) {
            dt = inputParams.tStp;
        }

        timeStepCount += 1;
        time += dt;

        V.divergence(divV);
        maxDivergence = divV.F.fieldMax();

        int iY = 0;
        localEnergy = 0.0;
        totalEnergy = 0.0;
        for (int iX = xLow; iX <= xTop; iX++) {
            for (int iZ = zLow; iZ <= zTop; iZ++) {
                localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*dVol;
            }
        }

        MPI_Allreduce(&localEnergy, &totalEnergy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);

        if (timeStepCount % inputParams.ioCnt == 0) {
            if (mesh.rankData.rank == 0) {
                std::cout << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                           std::setw(16) << std::setprecision(8) << sqrt(totalEnergy) << "\t" << maxDivergence << std::endl;

                ofFile << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                        std::setw(16) << std::setprecision(8) << sqrt(totalEnergy) << "\t" << maxDivergence << "\t" << dt << std::endl;
            }
        }

        if (inputParams.readProbes and std::abs(prTime - time) < 0.5*dt) {
            dataProbe->probeData(time);
            prTime += inputParams.prInt;
        }

        if (std::abs(fwTime - time) < 0.5*dt) {
            dataWriter.writeData(time);
            fwTime += inputParams.fwInt;
        }

        if (std::abs(inputParams.tMax - time) < 0.5*dt) {
            break;
        }
    }

    // CLOSE THE TIME-SERIES FILE
    ofFile.close();
}

void hydro_d2::computeTimeStep() {
    Hv = 0.0;

    // CALCULATE Hv FROM THE NON LINEAR TERMS AND HALF THE VISCOUS TERMS
    V.computeDiff(Hv);
    Hv *= inverseRe;

    // COMPUTE THE CONVECTIVE DERIVATIVE AND SUBTRACT IT FROM THE CALCULATED DIFFUSION TERMS OF RHS IN Hv
    V.computeNLin(V, Hv);

    Force.add_VForce(Hv);

    pressureGradient = 0.0;
    P.gradient(pressureGradient, V);

    // ADD PRESSURE GRADIENT TO NON-LINEAR TERMS AND MULTIPLY WITH TIME-STEP
    Hv -= pressureGradient;
    Hv *= dt;

    // ADD THE CALCULATED VALUES TO THE VELOCITY AT START OF TIME-STEP
    Hv += V;

    // SYNCHRONISE THE RHS OF TIME INTEGRATION STEP THUS OBTAINED ACROSS ALL PROCESSORS
    Hv.syncData();

    // CALCULATE V IMPLICITLY USING THE JACOBI ITERATIVE SOLVER
    solveVx();
    solveVz();

    // CALCULATE THE RHS FOR THE POISSON SOLVER FROM THE GUESSED VALUES OF VELOCITY IN V
    V.divergence(mgRHS);
    mgRHS *= 1.0/dt;

    // USING THE CALCULATED mgRHS, EVALUATE Pp USING MULTI-GRID METHOD
    mgSolver.mgSolve(Pp, mgRHS);

    // SYNCHRONISE THE PRESSURE CORRECTION ACROSS PROCESSORS
    Pp.syncData();

    // ADD THE PRESSURE CORRECTION CALCULATED FROM THE POISSON SOLVER TO P
    P += Pp;

    // CALCULATE FINAL VALUE OF V BY SUBTRACTING THE GRADIENT OF PRESSURE CORRECTION
    Pp.gradient(pressureGradient, V);
    pressureGradient *= dt;
    V -= pressureGradient;

    // IMPOSE BOUNDARY CONDITIONS ON V
    imposeUBCs();
    imposeWBCs();
}

void hydro_d2::solveVx() {
    int iterCount = 0;
    double maxError = 0.0;

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                guessedVelocity.Vx(iX, iY, iZ) = ((hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                                       hx2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                                 dt / ( hz2hx2 * 2.0 * inputParams.Re) + Hv.Vx.F(iX, iY, iZ))/
                             (1.0 + dt * ((hz2 * mesh.xix2Colloc(iX) +
                                                       hx2 * mesh.ztz2Staggr(iZ)))/(inputParams.Re * hz2hx2));
            }
        }

        V.Vx.F = guessedVelocity.Vx;

        imposeUBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                velocityLaplacian.Vx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - (
                                mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx2) +
                                mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz2)) *
                                            0.5 * dt * inverseRe;
            }
        }

        velocityLaplacian.Vx(V.Vx.fBulk) = abs(velocityLaplacian.Vx(V.Vx.fBulk) - Hv.Vx.F(V.Vx.fBulk));

        maxError = velocityLaplacian.vxMax();

        if (maxError < inputParams.tolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vx not converging. Aborting" << std::endl;
                exit(0);
            }
            break;
        }
    }
}

void hydro_d2::solveVz() {
    int iterCount = 0;
    double maxError = 0.0;

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                guessedVelocity.Vz(iX, iY, iZ) = ((hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                                       hx2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                                    dt / ( hz2hx2 * 2.0 * inputParams.Re) + Hv.Vz.F(iX, iY, iZ))/
                             (1.0 + dt * ((hz2 * mesh.xix2Staggr(iX) +
                                                       hx2 * mesh.ztz2Colloc(iZ)))/(inputParams.Re * hz2hx2));
            }
        }

        V.Vz.F = guessedVelocity.Vz;

        imposeWBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                velocityLaplacian.Vz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - (
                                mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx2) +
                                mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz2)) *
                                                0.5 * dt * inverseRe;
            }
        }

        velocityLaplacian.Vz(V.Vz.fBulk) = abs(velocityLaplacian.Vz(V.Vz.fBulk) - Hv.Vz.F(V.Vz.fBulk));

        maxError = velocityLaplacian.vzMax();

        if (maxError < inputParams.tolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vz not converging. Aborting" << std::endl;
                exit(0);
            }
            break;
        }
    }
}

void hydro_d2::setCoefficients() {
    hx = mesh.dXi;
    hz = mesh.dZt;

    hx2 = pow(mesh.dXi, 2.0);
    hz2 = pow(mesh.dZt, 2.0);

    hz2hx2 = pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose global and sub-domain boundary values on x-velocity component
 *
 *          First, inter-processor data transfer is performed by calling the \ref sfield#syncData "syncData" function of Vx
 *          Then the values of <B>Vx</B> on the 4 walls are imposed.
 *          The order of imposing boundary conditions is - left, right, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 ********************************************************************************************************************************************
 */
void hydro_d2::imposeUBCs() {
    V.Vx.syncData();

    // IMPOSE BC FOR Vx ALONG LEFT AND RIGHT WALLS
    if (not inputParams.xPer) {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1) {
            // Vx LIES ON EITHER SIDE OF THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
            if (mesh.rankData.xRank == 0) {
                V.Vx.F(V.Vx.fWalls(0)) = -V.Vx.F(V.Vx.xDim.shift(V.Vx.fWalls(0), 1));
            }
            // Vx LIES ON EITHER SIDE OF THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                V.Vx.F(V.Vx.fWalls(1)) = -V.Vx.F(V.Vx.xDim.shift(V.Vx.fWalls(1), -1));
            }
        // INFLOW AND OUTFLOW BCS
        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
            // Vx LIES ON EITHER SIDE OF THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
            if (mesh.rankData.xRank == 0) {
                // INFLOW BOUNDARY CONDITION AT INLET - IMPOSE VELOCITY OF 1.0
                V.Vx.F(V.Vx.fWalls(0)) = 2.0 - V.Vx.F(V.Vx.xDim.shift(V.Vx.fWalls(0), 1));
            }
            // Vx LIES ON EITHER SIDE OF THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                // OUTFLOW BOUNDARY CONDITION AT EXIT - NEUMANN BC WITH DERIVATIVE ALONG X SET TO 0
                V.Vx.F(V.Vx.fWalls(1)) = V.Vx.F(V.Vx.xDim.shift(V.Vx.fWalls(1), -1));
            }
        }
    } // FOR PERIODIC BCS, THE MPI DATA TRANSFER IS SUFFICIENT FOR COLLOCATED GRID POINTS AT LEFT AND RIGHT WALLS

    // IMPOSE BC FOR Vx ALONG TOP AND BOTTOM WALLS
    if (inputParams.zPer) {
        // PERIODIC BCS
        // Vx LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
        V.Vx.F(V.Vx.fWalls(4)) = 0.5*(V.Vx.F(V.Vx.zDim.shift(V.Vx.fWalls(4), -1)) + V.Vx.F(V.Vx.zDim.shift(V.Vx.fWalls(4), 1)));
        // Vx LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
        V.Vx.F(V.Vx.fWalls(5)) = 0.5*(V.Vx.F(V.Vx.zDim.shift(V.Vx.fWalls(5), -1)) + V.Vx.F(V.Vx.zDim.shift(V.Vx.fWalls(5), 1)));
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1) {
            // Vx LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
            V.Vx.F(V.Vx.fWalls(4)) = 0.0;
            // Vx LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
            V.Vx.F(V.Vx.fWalls(5)) = 1.0;
        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
            // Vx LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
            V.Vx.F(V.Vx.fWalls(4)) = 0.0;
            // Vx LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
            V.Vx.F(V.Vx.fWalls(5)) = 0.0;
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose global and sub-domain boundary values on z-velocity component
 *
 *          First, inter-processor data transfer is performed by calling the \ref sfield#syncData "syncData" function of Vz
 *          Then the values of <B>Vz</B> on the 4 walls are imposed.
 *          The order of imposing boundary conditions is - left, right, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 ********************************************************************************************************************************************
 */
void hydro_d2::imposeWBCs() {
    V.Vz.syncData();

    // IMPOSE BC FOR Vz ALONG LEFT AND RIGHT WALLS
    if (inputParams.xPer) {
        // PERIODIC BCS
        // Vz LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.xRank == 0) {
            V.Vz.F(V.Vz.fWalls(0)) = 0.5*(V.Vz.F(V.Vz.xDim.shift(V.Vz.fWalls(0), -1)) + V.Vz.F(V.Vz.xDim.shift(V.Vz.fWalls(0), 1)));
        }
        // Vz LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            V.Vz.F(V.Vz.fWalls(1)) = 0.5*(V.Vz.F(V.Vz.xDim.shift(V.Vz.fWalls(1), -1)) + V.Vz.F(V.Vz.xDim.shift(V.Vz.fWalls(1), 1)));
        }
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1) {
            // Vz LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
            if (mesh.rankData.xRank == 0) {
                V.Vz.F(V.Vz.fWalls(0)) = 0.0;
            }
            // Vz LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                V.Vz.F(V.Vz.fWalls(1)) = 0.0;
            }
        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
            // Vz LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
            if (mesh.rankData.xRank == 0) {
                // INFLOW BOUNDARY CONDITION AT INLET - IMPOSE VELOCITY OF 0.0
                V.Vz.F(V.Vz.fWalls(0)) = 0.0;
            }
            // Vz LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                // OUTFLOW BOUNDARY CONDITION AT EXIT - NEUMANN BC WITH DERIVATIVE ALONG X SET TO 0
                V.Vz.F(V.Vz.fWalls(1)) = V.Vz.F(V.Vz.zDim.shift(V.Vz.fWalls(1), -1));
            }
        }
    }

    // IMPOSE BC FOR Vz ALONG TOP AND BOTTOM WALLS
    if (inputParams.zPer) {
        // PERIODIC BCS
        V.Vz.F(V.Vz.fWalls(4)) = V.Vz.F(V.Vz.zDim.shift(V.Vz.fWalls(5), -1));
        V.Vz.F(V.Vz.fWalls(5)) = V.Vz.F(V.Vz.zDim.shift(V.Vz.fWalls(4), 1));
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1 or inputParams.probType == 3 or inputParams.probType == 4) {
            // Vz LIES ON EITHER SIDE OF THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS COLLOCATED ALONG Z
            V.Vz.F(V.Vz.fWalls(4)) = -V.Vz.F(V.Vz.zDim.shift(V.Vz.fWalls(4), 1));
            // Vz LIES ON EITHER SIDE OF THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS COLLOCATED ALONG Z
            V.Vz.F(V.Vz.fWalls(5)) = -V.Vz.F(V.Vz.zDim.shift(V.Vz.fWalls(5), -1));
        }
    }
}

double hydro_d2::testPeriodic() {
    int iY = 0;
    double xCoord = 0.0;
    double zCoord = 0.0;

    Hv = 0.0;
    V = 0.0;

    for (int i=V.Vx.F.lbound(0); i <= V.Vx.F.ubound(0); i++) {
        for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
            V.Vx.F(i, 0, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            Hv.Vx.F(i, 0, k) = V.Vx.F(i, 0, k);
        }
    }
    for (int i=V.Vz.F.lbound(0); i <= V.Vz.F.ubound(0); i++) {
        for (int k=V.Vz.F.lbound(2); k <= V.Vz.F.ubound(2); k++) {
            V.Vz.F(i, 0, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                 sin(2.0*M_PI*mesh.zColloc(k)/mesh.zLen);
            Hv.Vz.F(i, 0, k) = V.Vz.F(i, 0, k);
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS OF Vx IF DATA TRANSFER HAPPENS WITH NO HITCH
    // X-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += 1) {
            xCoord = mesh.xColloc(V.Vx.fCore.lbound(0)) - (mesh.xColloc(V.Vx.fCore.lbound(0) + iX) - mesh.xColloc(V.Vx.fCore.lbound(0)));
            Hv.Vx.F(V.Vx.fCore.lbound(0) - iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                         cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

            xCoord = mesh.xColloc(V.Vx.fCore.ubound(0)) + (mesh.xColloc(V.Vx.fCore.ubound(0)) - mesh.xColloc(V.Vx.fCore.ubound(0) - iX));
            Hv.Vx.F(V.Vx.fCore.ubound(0) + iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                         cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
        }
    }

    // X-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vx.fCore.lbound(0); iX <= V.Vx.fCore.ubound(0); iX += 1) {
            zCoord = mesh.zStaggr(V.Vx.fCore.lbound(2)) - (mesh.zStaggr(V.Vx.fCore.lbound(2) + iZ) - mesh.zStaggr(V.Vx.fCore.lbound(2)));
            Hv.Vx.F(iX, iY, V.Vx.fCore.lbound(2) - iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                         cos(2.0*M_PI*zCoord/mesh.zLen);

            zCoord = mesh.zStaggr(V.Vx.fCore.ubound(2)) + (mesh.zStaggr(V.Vx.fCore.ubound(2)) - mesh.zStaggr(V.Vx.fCore.ubound(2) - iZ));
            Hv.Vx.F(iX, iY, V.Vx.fCore.ubound(2) + iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                         cos(2.0*M_PI*zCoord/mesh.zLen);
        }
    }

    // Z-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iZ = V.Vz.fCore.lbound(2); iZ <= V.Vz.fCore.ubound(2); iZ += 1) {
            xCoord = mesh.xStaggr(V.Vz.fCore.lbound(0)) - (mesh.xStaggr(V.Vz.fCore.lbound(0) + iX) - mesh.xStaggr(V.Vz.fCore.lbound(0)));
            Hv.Vz.F(V.Vz.fCore.lbound(0) - iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                          sin(2.0*M_PI*mesh.zColloc(iZ)/mesh.zLen);

            xCoord = mesh.xStaggr(V.Vz.fCore.ubound(0)) + (mesh.xStaggr(V.Vz.fCore.ubound(0)) - mesh.xStaggr(V.Vz.fCore.ubound(0) - iX));
            Hv.Vz.F(V.Vz.fCore.ubound(0) + iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                          sin(2.0*M_PI*mesh.zColloc(iZ)/mesh.zLen);
        }
    }

    // Z-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vz.fCore.lbound(0); iX <= V.Vz.fCore.ubound(0); iX += 1) {
            zCoord = mesh.zColloc(V.Vz.fCore.lbound(2)) - (mesh.zColloc(V.Vz.fCore.lbound(2) + iZ) - mesh.zColloc(V.Vz.fCore.lbound(2)));
            Hv.Vz.F(iX, iY, V.Vz.fCore.lbound(2) - iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                          sin(2.0*M_PI*zCoord/mesh.zLen);

            zCoord = mesh.zColloc(V.Vz.fCore.ubound(2)) + (mesh.zColloc(V.Vz.fCore.ubound(2)) - mesh.zColloc(V.Vz.fCore.ubound(2) - iZ));
            Hv.Vz.F(iX, iY, V.Vz.fCore.ubound(2) + iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                          sin(2.0*M_PI*zCoord/mesh.zLen);
        }
    }

    imposeUBCs();
    imposeWBCs();

    V -= Hv;

    return std::max(blitz::max(fabs(V.Vx.F)), blitz::max(fabs(V.Vz.F)));
}

hydro_d2::~hydro_d2() { }
