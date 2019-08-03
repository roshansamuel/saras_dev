#include <iostream>
#include "parallel.h"
#include "scalar.h"
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
                std::cout << std::endl << "Solving NSE for heated bottom-plate problem" << std::endl;
            } else if (inputParams.probType == 6) {
                std::cout << std::endl << "Solving NSE for heated top-plate problem" << std::endl;
            } else {
                if (not inputParams.xPer) {
                    std::cout << std::endl << "Solving NSE for heated sidewall problem" << std::endl;
                } else {
                    std::cout << std::endl << "ERROR: X direction cannot be periodic for heated sidewall problem. ABORTING" << std::endl;

                    MPI_Finalize();
                    exit(0);
                }
            }
            std::cout << std::endl;
        }

        scalar *nseSolver;

        // CREATE NEW INSTANCE OF THE HYDRODYNAMICS SOLVER WITH SCALAR SOLVER
#ifdef PLANAR
        nseSolver = new scalar_d2(gridData, inputParams, mpi);
#else
        nseSolver = new scalar_d3(gridData, inputParams, mpi);
#endif

        nseSolver->solvePDE();

        delete nseSolver;

    }    
    
    else {
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
