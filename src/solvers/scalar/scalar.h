#ifndef SCALAR_H
#define SCALAR_H

#include <blitz/array.h>

#include "hydro.h"
#include <math.h>

class scalar: public hydro {
    public:
        sfield T;

        double nu, kappa; 

        scalar(const grid &mesh, const parser &solParam, parallel &mpiParam);

        virtual ~scalar() { };

    protected:
        boundary *tempBC;

        plainsf tmpRHS;

        plainsf guessedScalar;
        plainsf scalarLaplacian;

        virtual void solveT();
};

/**
 ********************************************************************************************************************************************
 *  \class scalar scalar.h "lib/scalar.h"
 *  \brief The base class scalar to solve the incompressible Navier-Stokes equations with energy equation
 *
 *  The class initializes and stores the velocity vector field and the pressure scalar field along with a few auxilliary
 *  fields to solve the PDE.
 *  In addition to the fields in hydro class, the scalar equation is also solved here.
 *  It solves the NSE using the \ref solvePDE function from within which the implicit Crank-Nicholson method is used
 *  to solve the PDE.
 ********************************************************************************************************************************************
 */

class scalar_d2: public scalar {
    public:
        scalar_d2(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        double testPeriodic();

        ~scalar_d2();

    private:
        multigrid_d2 mgSolver;

        void solveVx();
        void solveVz();

        void solveT();

        void imposeUBCs();
        void imposeWBCs();

        void computeTimeStep();
};

/**
 ********************************************************************************************************************************************
 *  \class scalar_d2 scalar.h "lib/scalar.h"
 *  \brief The derived class from the scalar base class to solve the incompressible NSE in 2D with energy equation
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Since the class is instantiated when solveing the NSE in 2D, the y-direction component of the grid is supressed.
 *  Consequently, the boundary conditions are imposed only on 4 sides of the domain.
 ********************************************************************************************************************************************
 */

class scalar_d3: public scalar {
    public:
        scalar_d3(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        double testPeriodic();

        ~scalar_d3();

    private:
        multigrid_d3 mgSolver;

#ifdef TIME_RUN
        double visc_time, nlin_time, intr_time, impl_time, prhs_time, pois_time;
#endif

        void solveVx();
        void solveVy();
        void solveVz();

        void solveT();

        void imposeUBCs();
        void imposeVBCs();
        void imposeWBCs();

        void computeTimeStep();
};

/**
 ********************************************************************************************************************************************
 *  \class scalar_d3 scalar.h "lib/scalar.h"
 *  \brief The derived class from the scalar base class to solve the incompressible NSE in 3D with energy equation
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Moreover, it imposes boundary conditions on all the three faces of the computational domain.
 ********************************************************************************************************************************************
 */

#endif
