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
/*! \file timestep.h
 *
 *  \brief Class declaration of timestep
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <blitz/array.h>

#include "force.h"
#include "poisson.h"

class timestep {
    public:
        timestep(const grid &mesh, const real &dt, vfield &V, sfield &P);

        virtual void timeAdvance(vfield &V, sfield &P);

    protected:
        const real &dt;

        const grid &mesh;

        real inverseRe;

        /** Plain scalar field into which the pressure correction is calculated and written by the Poisson solver */
        plainsf Pp;
        /** Plain scalar field into which the RHS for pressure Poisson equation is written and passed to the Poisson solver */
        plainsf mgRHS;

        /** Plain vector field into which the RHS of the Navier-Stokes equation is written and stored */
        plainvf nseRHS;
        /** Plain vector field which stores the pressure gradient term. */
        plainvf pressureGradient;
};

/**
 ********************************************************************************************************************************************
 *  \class timestep timestep.h "lib/timestep/timestep.h"
 *  \brief Contains all the global variables related to time-advancing the solution by one time-step
 *
 ********************************************************************************************************************************************
 */

class eulerCN_d2: public timestep {
    public:
        eulerCN_d2(const grid &mesh, const real &dt, vfield &V, sfield &P);

        void timeAdvance(vfield &V, sfield &P);

    private:
        /** Maximum number of iterations for the iterative solvers \ref hydro#solveVx, \ref hydro#solveVy and \ref hydro#solveVz */
        int maxIterations;

        real hx2, hz2, hz2hx2;

        /** Plain vector field which serves as a temporary array during iterative solution procedure for velocity terms. */
        plainvf guessedVelocity;

        /** Plain vector field into which the RHS of the implicit equation for velocities are calculated during iterative solving. */
        plainvf velocityLaplacian;

        multigrid_d2 mgSolver;

        void solveVx(vfield &V);
        void solveVz(vfield &V);

        void setCoefficients();
};

/**
 ********************************************************************************************************************************************
 *  \class eulerCN_d2 timestep.h "lib/timestep/timestep.h"
 *  \brief The derived class from timestep to advance the 2D solution using explicit Euler and implicit Crank-Nicholson methods
 *
 ********************************************************************************************************************************************
 */

class eulerCN_d3: public timestep {
    public:
        eulerCN_d3(const grid &mesh, const real &dt, vfield &V, sfield &P);

        void timeAdvance(vfield &V, sfield &P);

    private:
        /** Maximum number of iterations for the iterative solvers \ref hydro#solveVx, \ref hydro#solveVy and \ref hydro#solveVz */
        int maxIterations;

        real hx2, hy2, hz2;
        real hx2hy2, hz2hx2, hy2hz2, hx2hy2hz2;

#ifdef TIME_RUN
        real visc_time, nlin_time, intr_time, impl_time, prhs_time, pois_time;
#endif

        /** Plain vector field which serves as a temporary array during iterative solution procedure for velocity terms. */
        plainvf guessedVelocity;

        /** Plain vector field into which the RHS of the implicit equation for velocities are calculated during iterative solving. */
        plainvf velocityLaplacian;

        multigrid_d3 mgSolver;

        void solveVx(vfield &V);
        void solveVy(vfield &V);
        void solveVz(vfield &V);

        void setCoefficients();
};

/**
 ********************************************************************************************************************************************
 *  \class eulerCN_d3 timestep.h "lib/timestep/timestep.h"
 *  \brief The derived class from timestep to advance the 3D solution using explicit Euler and implicit Crank-Nicholson methods
 *
 ********************************************************************************************************************************************
 */

#endif
