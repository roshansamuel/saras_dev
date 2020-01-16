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
/*! \file hydro.h
 *
 *  \brief Class declaration of the hydro solver for both 2D and 3D cases.
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef HYDRO_H
#define HYDRO_H

#include <blitz/array.h>

#include "boundary.h"
#include "parallel.h"
#include "poisson.h"
#include "plainvf.h"
#include "tseries.h"
#include "writer.h"
#include "reader.h"
#include "probes.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"
#include "force.h"
#include "grid.h"

class hydro {
    public:
        vfield V;

        sfield P;

        force *vForcing;

        hydro(const grid &mesh, const parser &solParam, parallel &mpiParam);

        virtual void solvePDE();
        virtual double testPeriodic();

        virtual ~hydro() { };

    protected:
        int timeStepCount;
        int maxIterations;

        double time, dt;

        double hx, hy, hz;
        double hx2, hz2, hz2hx2;
        double hx2hy2, hy2hz2, hx2hy2hz2;

        const grid &mesh;
        const parser &inputParams;

        const double inverseRe;

        probes *dataProbe;
        boundary *uLft, *uRgt, *uFrn, *uBak, *uTop, *uBot;
        boundary *vLft, *vRgt, *vFrn, *vBak, *vTop, *vBot;
        boundary *wLft, *wRgt, *wFrn, *wBak, *wTop, *wBot;

        parallel &mpiData;

        plainsf Pp;
        plainsf mgRHS;

        plainvf nseRHS;
        plainvf velocityLaplacian;
        plainvf pressureGradient;
        plainvf guessedVelocity;

        void checkPeriodic();
        void setCoefficients();

        void initVBC();
        void imposeUBCs();
        void imposeVBCs();
        void imposeWBCs();

        void initVForcing();

        virtual void solveVx();
        virtual void solveVy();
        virtual void solveVz();

        virtual void computeTimeStep();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro hydro.h "lib/hydro.h"
 *  \brief The base class hydro to solve the incompressible Navier-Stokes equations
 *
 *  The class initializes and stores the velocity vector field and the pressure scalar field along with a few auxilliary
 *  fields to solve the PDE.
 *  It solves the NSE using the \ref solvePDE function from within which the implicit Crank-Nicholson method is used
 *  to solve the PDE.
 ********************************************************************************************************************************************
 */

class hydro_d2: public hydro {
    public:
        hydro_d2(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        double testPeriodic();

        ~hydro_d2();

    private:
        multigrid_d2 mgSolver;

        void solveVx();
        void solveVz();

        void computeTimeStep();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro_d2 hydro.h "lib/hydro.h"
 *  \brief The derived class from the hydro base class to solve the incompressible NSE in 2D
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Since the class is instantiated when solving the NSE in 2D, the y-direction component of the grid is supressed.
 *  Consequently, the boundary conditions are imposed only on 4 sides of the domain.
 ********************************************************************************************************************************************
 */

class hydro_d3: public hydro {
    public:
        hydro_d3(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        double testPeriodic();

        ~hydro_d3();

    private:
        multigrid_d3 mgSolver;

#ifdef TIME_RUN
        double visc_time, nlin_time, intr_time, impl_time, prhs_time, pois_time;
#endif

        void solveVx();
        void solveVy();
        void solveVz();

        void computeTimeStep();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro_d3 hydro.h "lib/hydro.h"
 *  \brief The derived class from the hydro base class to solve the incompressible NSE in 3D
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Moreover, it imposes boundary conditions on all the three faces of the computational domain.
 ********************************************************************************************************************************************
 */

#endif
