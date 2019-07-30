#include "thermal.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the base thermal class
 *
 *          The short base constructor of the thermal class merely assigns the const references to the grid and parser
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
thermal::thermal(const grid &mesh, const parser &solParam, parallel &mpiParam):
            hydro(mesh, solParam, mpiParam),
            guessedTemperature(mesh, "JAC_T", false),
            temperatureLaplacian(mesh, "LAP_T", false),
            Ht(mesh, "Ht", false),
            T(mesh, "T", true)
{
    // Below flags may be turned on for debugging/dignostic runs only
    bool viscSwitch = false;
    bool diffSwitch = false;

    if (inputParams.rbcType == 1) {
        nu = inputParams.Pr;
        kappa = 1.0;
        Fb = inputParams.Ra*inputParams.Pr;
    } else if (inputParams.rbcType == 2) {
        nu = sqrt(inputParams.Pr/inputParams.Ra);
        kappa = 1.0/sqrt(inputParams.Pr*inputParams.Ra);
        Fb = 1.0;
    } else if (inputParams.rbcType == 3) {
        nu = 1.0;
        kappa = 1.0/inputParams.Pr;
        Fb = inputParams.Ra;
    } else if (inputParams.rbcType == 4) {
        nu = sqrt(inputParams.Pr/inputParams.Ra);
        kappa = 1.0/sqrt(inputParams.Pr*inputParams.Ra);
        Fb = inputParams.Pr;
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
 * \brief   Function to solve the implicit equation for temperature
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
void thermal::solveT() { };
