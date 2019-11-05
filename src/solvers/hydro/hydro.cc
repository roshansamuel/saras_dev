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
 * \brief   Function to solve the implicit equation for x-velocity
 *
 *          The implicit equation for \f$ u_x' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
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
void hydro::solveVx() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for y-velocity
 *
 *          The implicit equation for \f$ u_y' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
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
void hydro::solveVy() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for z-velocity
 *
 *          The implicit equation for \f$ u_z' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
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
void hydro::solveVz() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the coefficients used for solving the implicit equations of U, V and W
 *
 *          The function assigns values to the variables \ref hx, \ref hy, etc.
 *          These coefficients are repeatedly used at many places in the Poisson solver for implicit calculation of velocities.
 ********************************************************************************************************************************************
 */
void hydro::setCoefficients() { };

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
 * \brief   Function to test whether periodic BC is being implemented properly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls imposeUBCs, imposeVBCs and imposeWBCs functions and checks if the correct values of the functions are imposed at boundaries
 ********************************************************************************************************************************************
 */
double hydro::testPeriodic() { return 0; };
