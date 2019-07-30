#include <iostream>
#include "parallel.h"
#include "thermal.h"
#include "parser.h"
#include "hydro.h"
#include "grid.h"

int main() {
    struct timeval runStart, runEnd;

    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    // ALL PROCESSES READ THE INPUT PARAMETERS
    parser inputParams;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputParams);

    grid gridData(inputParams, mpi);

    gettimeofday(&runStart, NULL);

    if (inputParams.probType <= 4) {
        if (mpi.rank == 0) {
            if (inputParams.probType == 1) {
                std::cout << std::endl << "Solving NSE for lid-driven cavity problem" << std::endl;
            } else if (inputParams.probType == 2) {
                std::cout << std::endl << "Solving NSE for Taylor-Green vortices" << std::endl;
            } else if (inputParams.probType == 3) {
                std::cout << std::endl << "Solving NSE for channel flow problem" << std::endl;
            } else if (inputParams.probType == 4) {
                std::cout << std::endl << "Solving NSE for square-duct flow problem" << std::endl;
            }
            std::cout << std::endl;
        }

        hydro *nseSolver;

        // CREATE NEW INSTANCE OF THE HYDRODYNAMICS SOLVER
#ifdef PLANAR
        nseSolver = new hydro_d2(gridData, inputParams, mpi);
#else
        nseSolver = new hydro_d3(gridData, inputParams, mpi);
#endif

        nseSolver->solvePDE();

        delete nseSolver;

    } else if (inputParams.probType <= 7) {
        if (mpi.rank == 0) {
            if (inputParams.probType == 5) {
                std::cout << std::endl << "Solving NSE for Rayleigh-Benard convection" << std::endl;
            } else if (inputParams.probType == 6) {
                std::cout << std::endl << "Solving NSE for Stably-stratified flows" << std::endl;
            } else {
                if (not inputParams.xPer) {
                    std::cout << std::endl << "Solving NSE for vertical convection" << std::endl;
                } else {
                    std::cout << std::endl << "ERROR: X direction cannot be periodic for vertical convection. ABORTING" << std::endl;

                    MPI_Finalize();
                    exit(0);
                }
            }
            std::cout << std::endl;
        }

        thermal *nseSolver;

        // CREATE NEW INSTANCE OF THE HYDRODYNAMICS SOLVER WITH TEMPERATURE SOLVER
#ifdef PLANAR
        nseSolver = new thermal_d2(gridData, inputParams, mpi);
#else
        nseSolver = new thermal_d3(gridData, inputParams, mpi);
#endif

        nseSolver->solvePDE();

        delete nseSolver;

    } else if (inputParams.probType == 8) {
#ifdef PLANAR
        if (mpi.rank == 0) {
            std::cout << std::endl << "ERROR: Rotating convection cannot be solved with PLANAR flag turned on. ABORTING" << std::endl;
        }

        MPI_Finalize();
        exit(0);
#endif
        if (mpi.rank == 0) {
            std::cout << std::endl << "Solving NSE for Rotating Rayleigh-Benard convection" << std::endl;
            std::cout << std::endl;
        }

        thermal *nseSolver;
        nseSolver = new thermal_d3(gridData, inputParams, mpi);

        nseSolver->solvePDE();

        delete nseSolver;

    } else {
        if (mpi.rank == 0) {
            std::cout << std::endl << "Invalid problem type. ABORTING" << std::endl;
        }

        MPI_Finalize();
        exit(0);
    }

    gettimeofday(&runEnd, NULL);
    double run_time = ((runEnd.tv_sec - runStart.tv_sec)*1000000u + runEnd.tv_usec - runStart.tv_usec)/1.e6;

    if (mpi.rank == 0) {
        std::cout << std::endl << "Simulation completed" << std::endl;
        std::cout << std::endl;
        std::cout << "Time taken by simulation: " << run_time << std::endl;
    }

    // FINALIZE AND CLEAN-UP
    MPI_Finalize();

    return 0;
}
