#ifndef THERMAL_H
#define THERMAL_H

#include <blitz/array.h>

#include "hydro.h"
#include <math.h>

class thermal: public hydro {
    protected:
        boundary *tempBC;

        sfield guessedTemperature;
        sfield temperatureLaplacian;

        sfield Ht;

        virtual void solveT();

    public:
        sfield T;

        double nu, kappa, Fb;

        thermal(const grid &mesh, const parser &solParam, parallel &mpiParam);

        virtual ~thermal() { };
};

/**
 ********************************************************************************************************************************************
 *  \class thermal thermal.h "lib/thermal.h"
 *  \brief The base class thermal to solve the incompressible Navier-Stokes equations with energy equation
 *
 *  The class initializes and stores the velocity vector field and the pressure scalar field along with a few auxilliary
 *  fields to solve the PDE.
 *  In addition to the fields in hydro class, the temperature equation is also solved here.
 *  It solves the NSE using the \ref solvePDE function from within which the implicit Crank-Nicholson method is used
 *  to solve the PDE.
 ********************************************************************************************************************************************
 */

class thermal_d2: public thermal {
    private:
        multigrid_d2 mgSolver;

        double hx2, hz2, hz2hx2;

        void solveVx();
        void solveVz();

        void solveT();

        void imposeUBCs();
        void imposeWBCs();

        void setCoefficients();

        void computeTimeStep();

    public:
        thermal_d2(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        double testPeriodic();

        ~thermal_d2();
};

/**
 ********************************************************************************************************************************************
 *  \class thermal_d2 thermal.h "lib/thermal.h"
 *  \brief The derived class from the thermal base class to solve the incompressible NSE in 2D with energy equation
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Since the class is instantiated when solveing the NSE in 2D, the y-direction component of the grid is supressed.
 *  Consequently, the boundary conditions are imposed only on 4 sides of the domain.
 ********************************************************************************************************************************************
 */

class thermal_d3: public thermal {
    private:
        multigrid_d3 mgSolver;

        double hx2hy2, hy2hz2, hz2hx2, hx2hy2hz2;

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

        void setCoefficients();

        void computeTimeStep();

    public:
        thermal_d3(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        double testPeriodic();

        ~thermal_d3();
};

/**
 ********************************************************************************************************************************************
 *  \class thermal_d3 thermal.h "lib/thermal.h"
 *  \brief The derived class from the thermal base class to solve the incompressible NSE in 3D with energy equation
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Moreover, it imposes boundary conditions on all the three faces of the computational domain.
 ********************************************************************************************************************************************
 */

#endif
