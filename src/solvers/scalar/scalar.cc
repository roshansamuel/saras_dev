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
/*! \file scalar.cc
 *
 *  \brief Definitions of common functions for both 2D and 3D runs of the solver class scalar - this class solves the basic Navier-Stokes equation.
 *  \sa scalar.h
 *  \author Roshan Samuel, Shashwat Bhattacharya
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "scalar.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the base scalar class
 *
 *          The short base constructor of the scalar class merely assigns the const references to the grid and parser
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
scalar::scalar(const grid &mesh, const parser &solParam, parallel &mpiParam):
            hydro(mesh, solParam, mpiParam),
            T(mesh, "T"),
            tmpRHS(mesh, T),
            guessedScalar(mesh, T),
            scalarLaplacian(mesh, T)
{
    // Below flags may be turned on for debugging/dignostic runs only
    bool viscSwitch = false;
    bool diffSwitch = false;

    if (inputParams.rbcType == 1) {
        nu = inputParams.Pr;
        kappa = 1.0;
    } else if (inputParams.rbcType == 2) {
        nu = sqrt(inputParams.Pr/inputParams.Ra);
        kappa = 1.0/sqrt(inputParams.Pr*inputParams.Ra);
    } else if (inputParams.rbcType == 3) {
        nu = 1.0;
        kappa = 1.0/inputParams.Pr;
    } else if (inputParams.rbcType == 4) {
        nu = sqrt(inputParams.Pr/inputParams.Ra);
        kappa = 1.0/sqrt(inputParams.Pr*inputParams.Ra);
    } else {
        if (mpiData.rank == 0) {
            std::cout << "ERROR: Invalid RBC non-dimensionalization type. Aborting" << std::endl;
        }
        exit(0);
    }

    // Additional option of turning off diffusion for debugging/diagnostics only
    if (viscSwitch) {
        nu = 0.0;
    }

    if (diffSwitch) {
        kappa = 0.0;
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for scalar field
 *
 *          The implicit equation for \f$ \theta' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref sfield#fieldMax "fieldMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void scalar::solveT() { };
