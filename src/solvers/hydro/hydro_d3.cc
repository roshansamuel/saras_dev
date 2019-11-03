#include <sys/time.h>
#include <ctime>

#include "hydro.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the hydro_d3 class derived from the base hydro class
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
hydro_d3::hydro_d3(const grid &mesh, const parser &solParam, parallel &mpiParam):
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
        readFields.push_back(V.Vy);
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
                for (int j=V.Vx.F.lbound(1); j <= V.Vx.F.ubound(1); j++) {
                    for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
                        V.Vx.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                            cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                            cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                    }
                }
            }

            // Y-VELOCITY
            for (int i=V.Vy.F.lbound(0); i <= V.Vy.F.ubound(0); i++) {
                for (int j=V.Vy.F.lbound(1); j <= V.Vy.F.ubound(1); j++) {
                    for (int k=V.Vy.F.lbound(2); k <= V.Vy.F.ubound(2); k++) {
                        V.Vy.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                            sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                            cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                    }
                }
            }

            // Z-VELOCITY
            V.Vz.F = 0.0;

        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
            // INITIALIZE PRESSURE TO 1.0 THROUGHOUT THE DOMAIN
            P = 1.0;

            // FOR CHANNEL AND DUCT FLOW, SET THE X-VELOCITY THROUGHOUT THE DOMAIN TO 1.0
            V.Vx.F = 1.0;
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

    // Front and rear walls
    if (not inputParams.yPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Y Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.yRank == 0)             mpiData.nearRanks(2) = MPI_PROC_NULL;
        if (mpiData.yRank == mpiData.npY-1) mpiData.nearRanks(3) = MPI_PROC_NULL;
    }

    // Inform user about BC along top and bottom walls
    if (not inputParams.zPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Z Direction" << std::endl;
            std::cout << std::endl;
        }
    }

    imposeUBCs();
    imposeVBCs();
    imposeWBCs();
}

void hydro_d3::solvePDE() {
#ifndef TIME_RUN
    int xLow, xTop;
    int yLow, yTop;
    int zLow, zTop;

    double dVol;
    double maxDivergence;
    double fwTime, prTime;
    double totalEnergy, localEnergy;

    // Scalar field to compute divergence
    sfield divV(mesh, "DIV_V", false);
#endif

#ifdef TIME_RUN
    visc_time = 0.0;
    nlin_time = 0.0;
    intr_time = 0.0;
    impl_time = 0.0;
    prhs_time = 0.0;
    pois_time = 0.0;

#else
    // Scalar fields to be written into HDF5 file is passed to writer class as a vector
    std::vector<field> writeFields;

    // Populate the vector with required scalar fields
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
    std::ofstream ofFile("output/TimeSeries.dat");

    // FILE WRITING TIME
    fwTime = time;

    // FIELD PROBING TIME
    prTime = time;
#endif

    timeStepCount = 0;

#ifndef TIME_RUN
    // PARAMETERS FOR COMPUTING TOTAL ENERGY IN THE DOMAIN
    // UPPER AND LOWER LIMITS WHEN COMPUTING ENERGY IN STAGGERED GRID
    xLow = P.F.fCore.lbound(0);        xTop = P.F.fCore.ubound(0);
    yLow = P.F.fCore.lbound(1);        yTop = P.F.fCore.ubound(1);
    zLow = P.F.fCore.lbound(2);        zTop = P.F.fCore.ubound(2);

    // STAGGERED GRIDS HAVE SHARED POINT ACROSS MPI-SUBDOMAINS - ACCORDINGLY DECREASE LIMITS
    if (mesh.rankData.xRank > 0) {
        xLow += 1;
    }

    if (mesh.rankData.yRank > 0) {
        yLow += 1;
    }

    // INFINITESIMAL VOLUME FOR INTEGRATING ENERGY OVER DOMAIN
    dVol = hx*hy*hz;

    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    V.divergence(divV);
    maxDivergence = divV.F.fieldMax();

    localEnergy = 0.0;
    totalEnergy = 0.0;
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iY = yLow; iY <= yTop; iY++) {
            for (int iZ = zLow; iZ <= zTop; iZ++) {
                localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vy.F(iX, iY-1, iZ) + V.Vy.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*dVol;
            }
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
#endif

    // TIME-INTEGRATION LOOP
    dt = inputParams.tStp;
    while (true) {
        // MAIN FUNCTION CALLED IN EACH LOOP TO UPDATE THE FIELDS AT EACH TIME-STEP
        computeTimeStep();
        V.computeTStp(dt);
        if (dt > inputParams.tStp) {
            dt = inputParams.tStp;
        }

        timeStepCount += 1;
        time += dt;

#ifndef TIME_RUN
        V.divergence(divV);
        maxDivergence = divV.F.fieldMax();

        localEnergy = 0.0;
        totalEnergy = 0.0;
        for (int iX = xLow; iX <= xTop; iX++) {
            for (int iY = yLow; iY <= yTop; iY++) {
                for (int iZ = zLow; iZ <= zTop; iZ++) {
                    localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                                    pow((V.Vy.F(iX, iY-1, iZ) + V.Vy.F(iX, iY, iZ))/2.0, 2.0) +
                                    pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*dVol;
                }
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

    // CLOSE THE TIME-SERIES FILE
#ifndef TIME_RUN
    ofFile.close();
#endif
}

void hydro_d3::computeTimeStep() {
#ifdef TIME_RUN
    struct timeval begin, end;
#endif

    Hv = 0.0;

    // CALCULATE Hv FROM THE NON LINEAR TERMS AND HALF THE VISCOUS TERMS
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    V.computeDiff(Hv);
    Hv *= inverseRe;
    gettimeofday(&end, NULL);
    visc_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    V.computeDiff(Hv);
    Hv *= inverseRe;
#endif

    // COMPUTE THE CONVECTIVE DERIVATIVE AND SUBTRACT IT FROM THE CALCULATED DIFFUSION TERMS OF RHS IN Hv
#ifndef TIME_RUN
    V.computeNLin(V, Hv);
#else
    gettimeofday(&begin, NULL);
    V.computeNLin(V, Hv);
    gettimeofday(&end, NULL);
    nlin_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

    gettimeofday(&begin, NULL);
#endif

    Force.add_VForce(Hv);

    pressureGradient = 0.0;
    P.gradient(pressureGradient, V);

    // ADD PRESSURE GRADIENT TO NON-LINEAR TERMS AND MULTIPLY WITH TIME-STEP
    Hv -= pressureGradient;
    Hv *= dt;

    // ADD THE CALCULATED VALUES TO THE VELOCITY AT START OF TIME-STEP
    Hv += V;

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    intr_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif

    // SYNCHRONISE THE RHS OF TIME INTEGRATION STEP THUS OBTAINED ACROSS ALL PROCESSORS
    Hv.syncData();

    // CALCULATE V IMPLICITLY USING THE JACOBI ITERATIVE SOLVER
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

    // CALCULATE THE RHS FOR THE POISSON SOLVER FROM THE GUESSED VALUES OF VELOCITY IN V
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    V.divergence(mgRHS);
    mgRHS *= 1.0/dt;
    gettimeofday(&end, NULL);
    prhs_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    V.divergence(mgRHS);
    mgRHS *= 1.0/dt;
#endif

    // USING THE CALCULATED mgRHS, EVALUATE Pp USING MULTI-GRID METHOD
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    mgSolver.mgSolve(Pp, mgRHS);
    gettimeofday(&end, NULL);
    pois_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    mgSolver.mgSolve(Pp, mgRHS);
#endif

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
    imposeVBCs();
    imposeWBCs();
}

void hydro_d3::solveVx() {
    int iterCount = 0;
    double maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vx(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                                           hz2hx2 * mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) + V.Vx.F(iX, iY-1, iZ)) +
                                                           hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                                     dt / ( hx2hy2hz2 * 2.0 * inputParams.Re) + Hv.Vx.F(iX, iY, iZ))/
                                 (1.0 + dt * ((hy2hz2 * mesh.xix2Colloc(iX) +
                                                           hz2hx2 * mesh.ety2Staggr(iY) +
                                                           hx2hy2 * mesh.ztz2Staggr(iZ)))/(inputParams.Re * hx2hy2hz2));
                }
            }
        }

        V.Vx.F = guessedVelocity.Vx;

        imposeUBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - (
                           mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx * hx) +
                           mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY-1, iZ)) / (hy * hy) +
                           mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz * hz)) *
                                           0.5 * dt * inverseRe;
                }
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

void hydro_d3::solveVy() {
    int iterCount = 0;
    double maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vy(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) +
                                                           hz2hx2 * mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) + V.Vy.F(iX, iY-1, iZ)) +
                                                           hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) + V.Vy.F(iX, iY, iZ-1))) *
                                     dt / ( hx2hy2hz2 * 2.0 * inputParams.Re) + Hv.Vy.F(iX, iY, iZ))/
                                 (1.0 + dt * ((hy2hz2 * mesh.xix2Staggr(iX) +
                                                           hz2hx2 * mesh.ety2Colloc(iY) +
                                                           hx2hy2 * mesh.ztz2Staggr(iZ)))/(inputParams.Re * hx2hy2hz2));
                }
            }
        }

        V.Vy.F = guessedVelocity.Vy;

        imposeVBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vy(iX, iY, iZ) = V.Vy.F(iX, iY, iZ) - (
                           mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) / (hx * hx) +
                           mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY-1, iZ)) / (hy * hy) +
                           mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY, iZ-1)) / (hz * hz)) *
                                           0.5 * dt * inverseRe;
                }
            }
        }

        velocityLaplacian.Vy(V.Vy.fBulk) = abs(velocityLaplacian.Vy(V.Vy.fBulk) - Hv.Vy.F(V.Vy.fBulk));

        maxError = velocityLaplacian.vyMax();

        if (maxError < inputParams.tolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vy not converging. Aborting" << std::endl;
                exit(0);
            }
            break;
        }
    }
}

void hydro_d3::solveVz() {
    int iterCount = 0;
    double maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vz(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                                           hz2hx2 * mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) + V.Vz.F(iX, iY-1, iZ)) +
                                                           hx2hy2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                                     dt / ( hx2hy2hz2 * 2.0 * inputParams.Re) + Hv.Vz.F(iX, iY, iZ))/
                                 (1.0 + dt * ((hy2hz2 * mesh.xix2Staggr(iX) +
                                                           hz2hx2 * mesh.ety2Staggr(iY) +
                                                           hx2hy2 * mesh.ztz2Colloc(iZ)))/(inputParams.Re * hx2hy2hz2));
                }
            }
        }

        V.Vz.F = guessedVelocity.Vz;

        imposeWBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - (
                           mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx * hx) +
                           mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY-1, iZ)) / (hy * hy) +
                           mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz * hz)) *
                                           0.5 * dt * inverseRe;
                }
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

void hydro_d3::setCoefficients() {
    hx = mesh.dXi;
    hy = mesh.dEt;
    hz = mesh.dZt;

    hx2hy2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0);
    hy2hz2 = pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
    hz2hx2 = pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);

    hx2hy2hz2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose global and sub-domain boundary values on x-velocity component
 *
 *          First, inter-processor data transfer is performed by calling the \ref sfield#syncData "syncData" function of Vx
 *          Then the values of <B>Vx</B> on the 6 walls are imposed.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 ********************************************************************************************************************************************
 */
void hydro_d3::imposeUBCs() {
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

    // IMPOSE BC FOR Vx ALONG FRONT AND BACK WALLS
    if (inputParams.yPer) {
        // PERIODIC BCS
        // Vx LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.yRank == 0) {
            V.Vx.F(V.Vx.fWalls(2)) = 0.5*(V.Vx.F(V.Vx.yDim.shift(V.Vx.fWalls(2), -1)) + V.Vx.F(V.Vx.yDim.shift(V.Vx.fWalls(2), 1)));
        }
        // Vx LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
            V.Vx.F(V.Vx.fWalls(3)) = 0.5*(V.Vx.F(V.Vx.yDim.shift(V.Vx.fWalls(3), -1)) + V.Vx.F(V.Vx.yDim.shift(V.Vx.fWalls(3), 1)));
        }
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1 or inputParams.probType == 4) {
            // Vx LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y
            if (mesh.rankData.yRank == 0) {
                V.Vx.F(V.Vx.fWalls(2)) = 0.0;
            }
            // Vx LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y
            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                V.Vx.F(V.Vx.fWalls(3)) = 0.0;
            }
        }
    }

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
 * \brief   Function to impose global and sub-domain boundary values on y-velocity component
 *
 *          First, inter-processor data transfer is performed by calling the \ref sfield#syncData "syncData" function of Vy
 *          Then the values of <B>Vy</B> on the 6 walls are imposed.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 ********************************************************************************************************************************************
 */
void hydro_d3::imposeVBCs() {
    V.Vy.syncData();

    // IMPOSE BC FOR Vy ALONG LEFT AND RIGHT WALLS
    if (inputParams.xPer) {
        // PERIODIC BCS
        // Vy LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.xRank == 0) {
            V.Vy.F(V.Vy.fWalls(0)) = 0.5*(V.Vy.F(V.Vy.xDim.shift(V.Vy.fWalls(0), -1)) + V.Vy.F(V.Vy.xDim.shift(V.Vy.fWalls(0), 1)));
        }
        // Vy LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            V.Vy.F(V.Vy.fWalls(1)) = 0.5*(V.Vy.F(V.Vy.xDim.shift(V.Vy.fWalls(1), -1)) + V.Vy.F(V.Vy.xDim.shift(V.Vy.fWalls(1), 1)));
        }
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1) {
            // Vy LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
            if (mesh.rankData.xRank == 0) {
                V.Vy.F(V.Vy.fWalls(0)) = 0.0;
            }
            // Vy LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                V.Vy.F(V.Vy.fWalls(1)) = 0.0;
            }
        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
            // Vy LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
            if (mesh.rankData.xRank == 0) {
                // INFLOW BOUNDARY CONDITION AT INLET - IMPOSE VELOCITY OF 0.0
                V.Vy.F(V.Vy.fWalls(0)) = 0.0;
            }
            // Vy LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                // OUTFLOW BOUNDARY CONDITION AT EXIT - NEUMANN BC WITH DERIVATIVE ALONG X SET TO 0
                V.Vy.F(V.Vy.fWalls(1)) = V.Vy.F(V.Vy.yDim.shift(V.Vy.fWalls(1), -1));
            }
        }
    }

    // IMPOSE BC FOR Vy ALONG FRONT AND BACK WALLS
    if (not inputParams.yPer) {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1) {
            // Vy LIES ON EITHER SIDE OF THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
            if (mesh.rankData.yRank == 0) {
                V.Vy.F(V.Vy.fWalls(2)) = -V.Vy.F(V.Vy.yDim.shift(V.Vy.fWalls(2), 1));
            }
            // Vy LIES ON EITHER SIDE OF THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                V.Vy.F(V.Vy.fWalls(3)) = -V.Vy.F(V.Vy.yDim.shift(V.Vy.fWalls(3), -1));
            }
        } else if (inputParams.probType == 4) {
            // Vy LIES ON EITHER SIDE OF THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
            if (mesh.rankData.yRank == 0) {
                V.Vy.F(V.Vy.fWalls(2)) = -V.Vy.F(V.Vy.yDim.shift(V.Vy.fWalls(2), 1));
            }
            // Vy LIES ON EITHER SIDE OF THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                V.Vy.F(V.Vy.fWalls(3)) = -V.Vy.F(V.Vy.yDim.shift(V.Vy.fWalls(3), -1));
            }
        }
    } // FOR PERIODIC BCS, THE MPI DATA TRANSFER IS SUFFICIENT FOR COLLOCATED GRID POINTS AT FRONT AND BACK WALLS

    // IMPOSE BC FOR Vy ALONG TOP AND BOTTOM WALLS
    if (inputParams.zPer) {
        // PERIODIC BCS
        // Vy LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
        V.Vy.F(V.Vy.fWalls(4)) = 0.5*(V.Vy.F(V.Vy.zDim.shift(V.Vy.fWalls(4), -1)) + V.Vy.F(V.Vy.zDim.shift(V.Vy.fWalls(4), 1)));
        // Vy LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
        V.Vy.F(V.Vy.fWalls(5)) = 0.5*(V.Vy.F(V.Vy.zDim.shift(V.Vy.fWalls(5), -1)) + V.Vy.F(V.Vy.zDim.shift(V.Vy.fWalls(5), 1)));
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1 or inputParams.probType == 3 or inputParams.probType == 4) {
            // Vy LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z
            V.Vy.F(V.Vy.fWalls(4)) = 0.0;
            // Vy LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z
            V.Vy.F(V.Vy.fWalls(5)) = 0.0;
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose global and sub-domain boundary values on z-velocity component
 *
 *          First, inter-processor data transfer is performed by calling the \ref sfield#syncData "syncData" function of Vz
 *          Then the values of <B>Vz</B> on the 6 walls are imposed.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 ********************************************************************************************************************************************
 */
void hydro_d3::imposeWBCs() {
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

    // IMPOSE BC FOR Vz ALONG FRONT AND BACK WALLS
    if (inputParams.yPer) {
        // PERIODIC BCS
        // Vz LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.yRank == 0) {
            V.Vz.F(V.Vz.fWalls(2)) = 0.5*(V.Vz.F(V.Vz.yDim.shift(V.Vz.fWalls(2), -1)) + V.Vz.F(V.Vz.yDim.shift(V.Vz.fWalls(2), 1)));
        }
        // Vz LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
            V.Vz.F(V.Vz.fWalls(3)) = 0.5*(V.Vz.F(V.Vz.yDim.shift(V.Vz.fWalls(3), -1)) + V.Vz.F(V.Vz.yDim.shift(V.Vz.fWalls(3), 1)));
        }
    } else {
        // NON PERIODIC BCS
        // NO-SLIP BCS
        if (inputParams.probType == 1 or inputParams.probType == 4) {
            // Vz LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y
            if (mesh.rankData.yRank == 0) {
                V.Vz.F(V.Vz.fWalls(2)) = 0.0;
            }
            // Vz LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y
            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
                V.Vz.F(V.Vz.fWalls(3)) = 0.0;
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

double hydro_d3::testPeriodic() {
    double xCoord = 0.0;
    double yCoord = 0.0;
    double zCoord = 0.0;

    Hv = 0.0;
    V = 0.0;

    for (int i=V.Vx.F.lbound(0); i <= V.Vx.F.ubound(0); i++) {
        for (int j=V.Vx.F.lbound(1); j <= V.Vx.F.ubound(1); j++) {
            for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
                V.Vx.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                     cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                     cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                Hv.Vx.F(i, j, k) = V.Vx.F(i, j, k);
            }
        }
    }

    for (int i=V.Vy.F.lbound(0); i <= V.Vy.F.ubound(0); i++) {
        for (int j=V.Vy.F.lbound(1); j <= V.Vy.F.ubound(1); j++) {
            for (int k=V.Vy.F.lbound(2); k <= V.Vy.F.ubound(2); k++) {
                V.Vy.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                      sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                      cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                Hv.Vy.F(i, j, k) = V.Vy.F(i, j, k);
            }
        }
    }

    V.Vz.F = 0.0;
    Hv.Vz.F = 0.0;

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    // X-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iY = V.Vx.fCore.lbound(1); iY <= V.Vx.fCore.ubound(1); iY += iX) {
            for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += iX) {
                xCoord = mesh.xColloc(V.Vx.fCore.lbound(0)) - (mesh.xColloc(V.Vx.fCore.lbound(0) + iX) - mesh.xColloc(V.Vx.fCore.lbound(0)));
                Hv.Vx.F(V.Vx.fCore.lbound(0) - iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                                  cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                                  cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xColloc(V.Vx.fCore.ubound(0)) + (mesh.xColloc(V.Vx.fCore.ubound(0)) - mesh.xColloc(V.Vx.fCore.ubound(0) - iX));
                Hv.Vx.F(V.Vx.fCore.ubound(0) + iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
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
                Hv.Vx.F(iX, V.Vx.fCore.lbound(1) - iY, iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                                  cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                                  cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yStaggr(V.Vx.fCore.ubound(1)) + (mesh.yStaggr(V.Vx.fCore.ubound(1)) - mesh.yStaggr(V.Vx.fCore.ubound(1) - iY));
                Hv.Vx.F(iX, V.Vx.fCore.ubound(1) + iY, iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
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
                Hv.Vx.F(iX, iY, V.Vx.fCore.lbound(2) - iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                                  cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                                  cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(V.Vx.fCore.ubound(2)) + (mesh.zStaggr(V.Vx.fCore.ubound(2)) - mesh.zStaggr(V.Vx.fCore.ubound(2) - iZ));
                Hv.Vx.F(iX, iY, V.Vx.fCore.ubound(2) + iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
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
                Hv.Vy.F(V.Vy.fCore.lbound(0) - iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                                   sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                   cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xStaggr(V.Vy.fCore.ubound(0)) + (mesh.xStaggr(V.Vy.fCore.ubound(0)) - mesh.xStaggr(V.Vy.fCore.ubound(0) - iX));
                Hv.Vy.F(V.Vy.fCore.ubound(0) + iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
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
                Hv.Vy.F(iX, V.Vy.fCore.lbound(1) - iY, iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                   sin(2.0*M_PI*yCoord/mesh.yLen)*
                                                                   cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yColloc(V.Vy.fCore.ubound(1)) + (mesh.yColloc(V.Vy.fCore.ubound(1)) - mesh.yColloc(V.Vy.fCore.ubound(1) - iY));
                Hv.Vy.F(iX, V.Vy.fCore.ubound(1) + iY, iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
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
                Hv.Vy.F(iX, iY, V.Vy.fCore.lbound(2) - iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                   sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                   cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(V.Vy.fCore.ubound(2)) + (mesh.zStaggr(V.Vy.fCore.ubound(2)) - mesh.zStaggr(V.Vy.fCore.ubound(2) - iZ));
                Hv.Vy.F(iX, iY, V.Vy.fCore.ubound(2) + iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                   sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                   cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    imposeUBCs();
    imposeVBCs();
    imposeWBCs();

    V -= Hv;

    return std::max(blitz::max(fabs(V.Vx.F)), blitz::max(fabs(V.Vy.F)));
}

hydro_d3::~hydro_d3() { }
