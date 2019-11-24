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
 *          Also, the maximum allowable number of iterations for the Jacobi iterative solver being used to solve for the
 *          velocities implicitly is set as \f$ N_{max} = N_x \times N_y \times N_z \f$, where \f$N_x\f$, \f$N_y\f$ and \f$N_z\f$
 *          are the number of grid points in the collocated grid at the local sub-domains along x, y and z directions
 *          respectively.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
hydro::hydro(const grid &mesh, const parser &solParam, parallel &mpiParam):
            V(mesh, "V"),
            P(mesh, "P"),
            Force(V, solParam, mpiParam),
            mesh(mesh),
            inputParams(solParam),
            inverseRe(1.0/inputParams.Re),
            mpiData(mpiParam),
            Pp(mesh, P),
            mgRHS(mesh, P),
            nseRHS(mesh, V),
            velocityLaplacian(mesh, V),
            pressureGradient(mesh, V),
            guessedVelocity(mesh, V)
{
    maxIterations = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);
}

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

        if (mpiData.xRank == 0)             mpiData.nearRanks(0) = MPI_PROC_NULL;
        if (mpiData.xRank == mpiData.npX-1) mpiData.nearRanks(1) = MPI_PROC_NULL;
    }

    // Front and rear walls
#ifdef PLANAR
    // Front and rear walls are by default non-periodic for 2D simulations
    if (mpiData.yRank == 0)             mpiData.nearRanks(2) = MPI_PROC_NULL;
    if (mpiData.yRank == mpiData.npY-1) mpiData.nearRanks(3) = MPI_PROC_NULL;

#else
    if (not inputParams.yPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Y Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.yRank == 0)             mpiData.nearRanks(2) = MPI_PROC_NULL;
        if (mpiData.yRank == mpiData.npY-1) mpiData.nearRanks(3) = MPI_PROC_NULL;
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
 * \brief   Function to solve the implicit equation for x-velocity
 *
 *          The implicit equation for \f$ u_x' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref plainvf#vxMax "vxMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void hydro::solveVx() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for y-velocity
 *
 *          The implicit equation for \f$ u_y' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref plainvf#vyMax "vyMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void hydro::solveVy() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for z-velocity
 *
 *          The implicit equation for \f$ u_z' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref plainvf#vzMax "vzMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void hydro::solveVz() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the coefficients used for solving the implicit equations of U, V and W
 *
 *          The function assigns values to the variables \ref hx, \ref hy, etc.
 *          These coefficients are repeatedly used at many places in the Poisson solver for implicit calculation of velocities.
 ********************************************************************************************************************************************
 */
void hydro::setCoefficients() {
    hx = mesh.dXi;
    hz = mesh.dZt;

    hz2hx2 = pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);

#ifdef PLANAR
    hx2 = pow(mesh.dXi, 2.0);
    hz2 = pow(mesh.dZt, 2.0);

#else
    hy = mesh.dEt;

    hx2hy2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0);
    hy2hz2 = pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);

    hx2hy2hz2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
#endif
};

/**
 ********************************************************************************************************************************************
 * \brief   The subroutine to solve the NS equations using the implicit Crank-Nicholson method
 *
 *          This function uses the values of velocity vector field and pressure scalar field, along with a specifed time-step
 *          to update the values of both fields by one time-step.
 *          Hence this function has to be repeatedly called in a loop from within the \ref solvePDE function to solve the equations.
 ********************************************************************************************************************************************
 */
void hydro::computeTimeStep() { };

/**
 ********************************************************************************************************************************************
 * \brief   The core publicly accessible function of the \ref hydro class to solve the Navier-Stokes equations
 *
 *          The NSE are integrated in time from within this function by calling \ref computeTimeStep in a loop.
 *          The function keeps track of the non-dimensional time with \ref time and number of iterations with \ref iterCount.
 *          Both these values are continuously incremented from within the loop, and finally, when \ref time has reached the
 *          user-ser value in \ref parser#tMax "tMax", the time-integration loop is broken and the program exits.
 ********************************************************************************************************************************************
 */
void hydro::solvePDE() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the boundary conditions for velocity
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::initVBC() {
    // NO-SLIP BCS
    if (inputParams.probType == 1 or inputParams.probType > 4) {
        uLft = new dirichletFC(mesh, V.Vx, 0, 0.0);
        uRgt = new dirichletFC(mesh, V.Vx, 1, 0.0);
    // INFLOW AND OUTFLOW BCS
    } else if (inputParams.probType == 3 or inputParams.probType == 4) {
        uLft = new dirichletFC(mesh, V.Vx, 0, 1.0);
        uRgt = new neumannFC(mesh, V.Vx, 1, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    if (inputParams.probType == 1 or inputParams.probType >= 4) {
        uFrn = new dirichletCC(mesh, V.Vx, 2, 0.0);
        uBak = new dirichletCC(mesh, V.Vx, 3, 0.0);
    }
#endif

    // NO-SLIP BCS FOR LDC
    if (inputParams.probType == 1) {
        uBot = new dirichletCC(mesh, V.Vx, 4, 0.0);
        uTop = new dirichletCC(mesh, V.Vx, 5, 1.0);
    // NO-SLIP BCS
    } else if (inputParams.probType >= 3) {
        uBot = new dirichletCC(mesh, V.Vx, 4, 0.0);
        uTop = new dirichletCC(mesh, V.Vx, 5, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    if (inputParams.probType == 1 or inputParams.probType > 4) {
        vLft = new dirichletCC(mesh, V.Vy, 0, 0.0);
        vRgt = new dirichletCC(mesh, V.Vy, 1, 0.0);
    // INFLOW AND OUTFLOW BCS
    } else if (inputParams.probType == 3 or inputParams.probType == 4) {
        vLft = new dirichletCC(mesh, V.Vy, 0, 0.0);
        vRgt = new neumannCC(mesh, V.Vy, 1, 0.0);
    }

    // NO-SLIP BCS
    if (inputParams.probType == 1 or inputParams.probType >= 4) {
        vFrn = new dirichletFC(mesh, V.Vy, 2, 0.0);
        vBak = new dirichletFC(mesh, V.Vy, 3, 0.0);
    }

    if (inputParams.probType == 1 or inputParams.probType >= 3) {
        vBot = new dirichletCC(mesh, V.Vy, 4, 0.0);
        vTop = new dirichletCC(mesh, V.Vy, 5, 0.0);
    }
#endif

    // NO-SLIP BCS
    if (inputParams.probType == 1 or inputParams.probType > 4) {
        wLft = new dirichletCC(mesh, V.Vz, 0, 0.0);
        wRgt = new dirichletCC(mesh, V.Vz, 1, 0.0);
    // INFLOW AND OUTFLOW BCS
    } else if (inputParams.probType == 3 or inputParams.probType == 4) {
        wLft = new dirichletCC(mesh, V.Vz, 0, 0.0);
        wRgt = new neumannCC(mesh, V.Vz, 1, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    if (inputParams.probType == 1 or inputParams.probType >= 4) {
        wFrn = new dirichletCC(mesh, V.Vz, 2, 0.0);
        wBak = new dirichletCC(mesh, V.Vz, 3, 0.0);
    }
#endif

    if (inputParams.probType == 1 or inputParams.probType >= 3) {
        wBot = new dirichletFC(mesh, V.Vz, 4, 0.0);
        wTop = new dirichletFC(mesh, V.Vz, 5, 0.0);
    }
};

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for temperature
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::imposeUBCs() {
    V.Vx.syncData();

    if (not inputParams.xPer) {
        uLft->imposeBC();
        uRgt->imposeBC();
    }
#ifndef PLANAR
    if (not inputParams.yPer) {
        uFrn->imposeBC();
        uBak->imposeBC();
    }
#endif
    uTop->imposeBC();
    uBot->imposeBC();
};

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for temperature
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::imposeVBCs() {
    V.Vy.syncData();

    if (not inputParams.xPer) {
        vLft->imposeBC();
        vRgt->imposeBC();
    }
#ifndef PLANAR
    if (not inputParams.yPer) {
        vFrn->imposeBC();
        vBak->imposeBC();
    }
#endif
    vTop->imposeBC();
    vBot->imposeBC();
};

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for temperature
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::imposeWBCs() {
    V.Vz.syncData();

    if (not inputParams.xPer) {
        wLft->imposeBC();
        wRgt->imposeBC();
    }
#ifndef PLANAR
    if (not inputParams.yPer) {
        wFrn->imposeBC();
        wBak->imposeBC();
    }
#endif
    wTop->imposeBC();
    wBot->imposeBC();
};

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
//void hydro::imposeUBCs() {
//    V.Vx.syncData();
//
//    // IMPOSE BC FOR Vx ALONG LEFT AND RIGHT WALLS
//    if (not inputParams.xPer) {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType > 4) {
//            // Vx LIES ON EITHER SIDE OF THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
//            if (mesh.rankData.xRank == 0) {
//                V.Vx.F(V.Vx.fWalls(0)) = -V.Vx.F(V.Vx.shift(0, V.Vx.fWalls(0), 1));
//            }
//            // Vx LIES ON EITHER SIDE OF THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
//            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//                V.Vx.F(V.Vx.fWalls(1)) = -V.Vx.F(V.Vx.shift(0, V.Vx.fWalls(1), -1));
//            }
//        // INFLOW AND OUTFLOW BCS
//        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
//            // Vx LIES ON EITHER SIDE OF THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
//            if (mesh.rankData.xRank == 0) {
//                // INFLOW BOUNDARY CONDITION AT INLET - IMPOSE VELOCITY OF 1.0
//                V.Vx.F(V.Vx.fWalls(0)) = 2.0 - V.Vx.F(V.Vx.shift(0, V.Vx.fWalls(0), 1));
//            }
//            // Vx LIES ON EITHER SIDE OF THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS COLLOCATED ALONG X
//            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//                // OUTFLOW BOUNDARY CONDITION AT EXIT - NEUMANN BC WITH DERIVATIVE ALONG X SET TO 0
//                V.Vx.F(V.Vx.fWalls(1)) = V.Vx.F(V.Vx.shift(0, V.Vx.fWalls(1), -1));
//            }
//        }
//    } // FOR PERIODIC BCS, THE MPI DATA TRANSFER IS SUFFICIENT FOR COLLOCATED GRID POINTS AT LEFT AND RIGHT WALLS
//
//#ifndef PLANAR
//    // IMPOSE BC FOR Vx ALONG FRONT AND BACK WALLS
//    if (inputParams.yPer) {
//        // PERIODIC BCS
//        // Vx LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.yRank == 0) {
//            V.Vx.F(V.Vx.fWalls(2)) = 0.5*(V.Vx.F(V.Vx.shift(1, V.Vx.fWalls(2), -1)) + V.Vx.F(V.Vx.shift(1, V.Vx.fWalls(2), 1)));
//        }
//        // Vx LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
//            V.Vx.F(V.Vx.fWalls(3)) = 0.5*(V.Vx.F(V.Vx.shift(1, V.Vx.fWalls(3), -1)) + V.Vx.F(V.Vx.shift(1, V.Vx.fWalls(3), 1)));
//        }
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType >= 4) {
//            // Vx LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y
//            if (mesh.rankData.yRank == 0) {
//                V.Vx.F(V.Vx.fWalls(2)) = 0.0;
//            }
//            // Vx LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Y
//            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
//                V.Vx.F(V.Vx.fWalls(3)) = 0.0;
//            }
//        }
//    }
//#endif
//
//    // IMPOSE BC FOR Vx ALONG TOP AND BOTTOM WALLS
//    if (inputParams.zPer) {
//        // PERIODIC BCS
//        // Vx LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
//        V.Vx.F(V.Vx.fWalls(4)) = 0.5*(V.Vx.F(V.Vx.shift(2, V.Vx.fWalls(4), -1)) + V.Vx.F(V.Vx.shift(2, V.Vx.fWalls(4), 1)));
//        // Vx LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
//        V.Vx.F(V.Vx.fWalls(5)) = 0.5*(V.Vx.F(V.Vx.shift(2, V.Vx.fWalls(5), -1)) + V.Vx.F(V.Vx.shift(2, V.Vx.fWalls(5), 1)));
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1) {
//            // Vx LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
//            V.Vx.F(V.Vx.fWalls(4)) = 0.0;
//            // Vx LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
//            V.Vx.F(V.Vx.fWalls(5)) = 1.0;
//        } else if (inputParams.probType >= 3) {
//            // Vx LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
//            V.Vx.F(V.Vx.fWalls(4)) = 0.0;
//            // Vx LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vx IS STAGGERED ALONG Z
//            V.Vx.F(V.Vx.fWalls(5)) = 0.0;
//        }
//    }
//};

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
//void hydro::imposeVBCs() {
//    V.Vy.syncData();
//
//    // IMPOSE BC FOR Vy ALONG LEFT AND RIGHT WALLS
//    if (inputParams.xPer) {
//        // PERIODIC BCS
//        // Vy LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.xRank == 0) {
//            V.Vy.F(V.Vy.fWalls(0)) = 0.5*(V.Vy.F(V.Vy.shift(0, V.Vy.fWalls(0), -1)) + V.Vy.F(V.Vy.shift(0, V.Vy.fWalls(0), 1)));
//        }
//        // Vy LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//            V.Vy.F(V.Vy.fWalls(1)) = 0.5*(V.Vy.F(V.Vy.shift(0, V.Vy.fWalls(1), -1)) + V.Vy.F(V.Vy.shift(0, V.Vy.fWalls(1), 1)));
//        }
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType > 4) {
//            // Vy LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == 0) {
//                V.Vy.F(V.Vy.fWalls(0)) = 0.0;
//            }
//            // Vy LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//                V.Vy.F(V.Vy.fWalls(1)) = 0.0;
//            }
//        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
//            // Vy LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == 0) {
//                // INFLOW BOUNDARY CONDITION AT INLET - IMPOSE VELOCITY OF 0.0
//                V.Vy.F(V.Vy.fWalls(0)) = 0.0;
//            }
//            // Vy LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//                // OUTFLOW BOUNDARY CONDITION AT EXIT - NEUMANN BC WITH DERIVATIVE ALONG X SET TO 0
//                V.Vy.F(V.Vy.fWalls(1)) = V.Vy.F(V.Vy.shift(1, V.Vy.fWalls(1), -1));
//            }
//        }
//    }
//
//    // IMPOSE BC FOR Vy ALONG FRONT AND BACK WALLS
//    if (not inputParams.yPer) {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType > 4) {
//            // Vy LIES ON EITHER SIDE OF THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
//            if (mesh.rankData.yRank == 0) {
//                V.Vy.F(V.Vy.fWalls(2)) = -V.Vy.F(V.Vy.shift(1, V.Vy.fWalls(2), 1));
//            }
//            // Vy LIES ON EITHER SIDE OF THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
//            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
//                V.Vy.F(V.Vy.fWalls(3)) = -V.Vy.F(V.Vy.shift(1, V.Vy.fWalls(3), -1));
//            }
//        } else if (inputParams.probType == 4) {
//            // Vy LIES ON EITHER SIDE OF THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
//            if (mesh.rankData.yRank == 0) {
//                V.Vy.F(V.Vy.fWalls(2)) = -V.Vy.F(V.Vy.shift(1, V.Vy.fWalls(2), 1));
//            }
//            // Vy LIES ON EITHER SIDE OF THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS COLLOCATED ALONG Y
//            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
//                V.Vy.F(V.Vy.fWalls(3)) = -V.Vy.F(V.Vy.shift(1, V.Vy.fWalls(3), -1));
//            }
//        }
//    } // FOR PERIODIC BCS, THE MPI DATA TRANSFER IS SUFFICIENT FOR COLLOCATED GRID POINTS AT FRONT AND BACK WALLS
//
//    // IMPOSE BC FOR Vy ALONG TOP AND BOTTOM WALLS
//    if (inputParams.zPer) {
//        // PERIODIC BCS
//        // Vy LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
//        V.Vy.F(V.Vy.fWalls(4)) = 0.5*(V.Vy.F(V.Vy.shift(2, V.Vy.fWalls(4), -1)) + V.Vy.F(V.Vy.shift(2, V.Vy.fWalls(4), 1)));
//        // Vy LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
//        V.Vy.F(V.Vy.fWalls(5)) = 0.5*(V.Vy.F(V.Vy.shift(2, V.Vy.fWalls(5), -1)) + V.Vy.F(V.Vy.shift(2, V.Vy.fWalls(5), 1)));
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType >= 3) {
//            // Vy LIES ON THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z
//            V.Vy.F(V.Vy.fWalls(4)) = 0.0;
//            // Vy LIES ON THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vy IS STAGGERED ALONG Z
//            V.Vy.F(V.Vy.fWalls(5)) = 0.0;
//        }
//    }
//};

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
//void hydro::imposeWBCs() {
//    V.Vz.syncData();
//
//    // IMPOSE BC FOR Vz ALONG LEFT AND RIGHT WALLS
//    if (inputParams.xPer) {
//        // PERIODIC BCS
//        // Vz LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.xRank == 0) {
//            V.Vz.F(V.Vz.fWalls(0)) = 0.5*(V.Vz.F(V.Vz.shift(0, V.Vz.fWalls(0), -1)) + V.Vz.F(V.Vz.shift(0, V.Vz.fWalls(0), 1)));
//        }
//        // Vz LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Z - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//            V.Vz.F(V.Vz.fWalls(1)) = 0.5*(V.Vz.F(V.Vz.shift(0, V.Vz.fWalls(1), -1)) + V.Vz.F(V.Vz.shift(0, V.Vz.fWalls(1), 1)));
//        }
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType > 4) {
//            // Vz LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == 0) {
//                V.Vz.F(V.Vz.fWalls(0)) = 0.0;
//            }
//            // Vz LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//                V.Vz.F(V.Vz.fWalls(1)) = 0.0;
//            }
//        } else if (inputParams.probType == 3 or inputParams.probType == 4) {
//            // Vz LIES ON THE LEFT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == 0) {
//                // INFLOW BOUNDARY CONDITION AT INLET - IMPOSE VELOCITY OF 0.0
//                V.Vz.F(V.Vz.fWalls(0)) = 0.0;
//            }
//            // Vz LIES ON THE RIGHT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG X
//            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
//                // OUTFLOW BOUNDARY CONDITION AT EXIT - NEUMANN BC WITH DERIVATIVE ALONG X SET TO 0
//                V.Vz.F(V.Vz.fWalls(1)) = V.Vz.F(V.Vz.shift(2, V.Vz.fWalls(1), -1));
//            }
//        }
//    }
//
//#ifndef PLANAR
//    // IMPOSE BC FOR Vz ALONG FRONT AND BACK WALLS
//    if (inputParams.yPer) {
//        // PERIODIC BCS
//        // Vz LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.yRank == 0) {
//            V.Vz.F(V.Vz.fWalls(2)) = 0.5*(V.Vz.F(V.Vz.shift(1, V.Vz.fWalls(2), -1)) + V.Vz.F(V.Vz.shift(1, V.Vz.fWalls(2), 1)));
//        }
//        // Vz LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y - PERIODIC BC IS IMPOSED BY AVERAGING
//        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
//            V.Vz.F(V.Vz.fWalls(3)) = 0.5*(V.Vz.F(V.Vz.shift(1, V.Vz.fWalls(3), -1)) + V.Vz.F(V.Vz.shift(1, V.Vz.fWalls(3), 1)));
//        }
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType >= 4) {
//            // Vz LIES ON THE FRONT WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y
//            if (mesh.rankData.yRank == 0) {
//                V.Vz.F(V.Vz.fWalls(2)) = 0.0;
//            }
//            // Vz LIES ON THE BACK WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS STAGGERED ALONG Y
//            if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
//                V.Vz.F(V.Vz.fWalls(3)) = 0.0;
//            }
//        }
//    }
//#endif
//
//    // IMPOSE BC FOR Vz ALONG TOP AND BOTTOM WALLS
//    if (inputParams.zPer) {
//        // PERIODIC BCS
//        V.Vz.F(V.Vz.fWalls(4)) = V.Vz.F(V.Vz.shift(2, V.Vz.fWalls(5), -1));
//        V.Vz.F(V.Vz.fWalls(5)) = V.Vz.F(V.Vz.shift(2, V.Vz.fWalls(4), 1));
//    } else {
//        // NON PERIODIC BCS
//        // NO-SLIP BCS
//        if (inputParams.probType == 1 or inputParams.probType >= 3) {
//            // Vz LIES ON EITHER SIDE OF THE BOTTOM WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS COLLOCATED ALONG Z
//            V.Vz.F(V.Vz.fWalls(4)) = -V.Vz.F(V.Vz.shift(2, V.Vz.fWalls(4), 1));
//            // Vz LIES ON EITHER SIDE OF THE TOP WALL AS THE WALL IS ON STAGGERED POINT AND Vz IS COLLOCATED ALONG Z
//            V.Vz.F(V.Vz.fWalls(5)) = -V.Vz.F(V.Vz.shift(2, V.Vz.fWalls(5), -1));
//        }
//    }
//};

/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether periodic BC is being implemented properly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls imposeUBCs, imposeVBCs and imposeWBCs functions and checks if the correct values of the functions are imposed at boundaries
 ********************************************************************************************************************************************
 */
double hydro::testPeriodic() { return 0; };
