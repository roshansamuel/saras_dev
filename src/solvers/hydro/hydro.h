#ifndef HYDRO_H
#define HYDRO_H

#include <blitz/array.h>

#include "boundary.h"
#include "parallel.h"
#include "poisson.h"
#include "plainvf.h"
#include "writer.h"
#include "reader.h"
#include "probes.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"
#include "grid.h"
#include "force.h"

class hydro {
    public:
        vfield V;

        sfield P;

        force Force;

        hydro(const grid &mesh, const parser &solParam, parallel &mpiParam);

        virtual void solvePDE();
        virtual double testPeriodic();

        virtual ~hydro() { };

    protected:
        int timeStepCount;
        int maxIterations;

        double time, dt;

        double hx, hy, hz;

        const grid &mesh;
        const parser &inputParams;

        const double inverseRe;

        probes *dataProbe;

        parallel &mpiData;

        plainsf Pp;
        plainsf mgRHS;

        plainvf nseRHS;
        plainvf velocityLaplacian;
        plainvf pressureGradient;
        plainvf guessedVelocity;

        virtual void solveVx();
        virtual void solveVy();
        virtual void solveVz();

        virtual void setCoefficients();

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

        double hx2, hz2, hz2hx2;

        void solveVx();
        void solveVz();

        void imposeUBCs();
        void imposeWBCs();

        void setCoefficients();

        void computeTimeStep();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro_d2 hydro.h "lib/hydro.h"
 *  \brief The derived class from the hydro base class to solve the incompressible NSE in 2D
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Since the class is instantiated when solveing the NSE in 2D, the y-direction component of the grid is supressed.
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

        double hx2hy2, hy2hz2, hz2hx2, hx2hy2hz2;

#ifdef TIME_RUN
        double visc_time, nlin_time, intr_time, impl_time, prhs_time, pois_time;
#endif

        void solveVx();
        void solveVy();
        void solveVz();

        void imposeUBCs();
        void imposeVBCs();
        void imposeWBCs();

        void setCoefficients();

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
