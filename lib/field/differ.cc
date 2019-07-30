#include "differ.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the differ class
 *
 *          The empty constructor of the differ class only serves to assign values to the two const parameters of the differ class,
 *          namely <B>dim</B> and <B>h</B>.
 *
 * \param   dim is an integer value that specifies along which direction the differ class' operations will act.
 * \param   h is a double precision value that specifies the uniform grid spacing to be used in finite differencing stencil
 ********************************************************************************************************************************************
 */
differ::differ(int dim, double h): dim(dim), h(h) { }

/**
 ********************************************************************************************************************************************
 * \brief   Function to shift a blitz RectDomain object a specified number of steps
 *
 *          The RectDomain objects offer a view of the blitz arrays on which the differ class operates.
 *          These objects are shifted along the dimension specified in the constructor, by <B>dim</B>, through a number of steps,
 *          to offer offset views for quick computation of finite differences.
 *
 * \param   core is the input RectDomain object which is to be shifted to get the new view
 * \param   steps is the integer value by which the input view must be offset along the dimension specified by <B>dim</B>
 *
 * \return  A RectDomain object that specifies the new offset view of the data
 ********************************************************************************************************************************************
 */
blitz::RectDomain<3> differ::shift(blitz::RectDomain<3> core, int steps) {
    core.lbound()(dim) += steps;
    core.ubound()(dim) += steps;

    return core;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the first derivative by finite-differencing stencil
 *
 *          The first derivative of the data specified through input to the function is computed using the second order,
 *          central finite-differencing stencil and written into an output array passed to the function.
 *
 * \param   inputMat is the 3D blitz array of which the derivatives have to be computed.
 * \param   outputMat is the 3D blitz array into which the calculated values of derivative are written
 * \param   core is the blitz RectDomain object which offers a view of the core of the matrix
 ********************************************************************************************************************************************
 */
void differ::D1D(blitz::Array<double, 3> inputMat, blitz::Array<double, 3> outputMat, blitz::RectDomain<3> core) {
    outputMat(core) = ((inputMat(shift(core, 1)) - inputMat(shift(core, -1)))/(2.0*h));
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the second derivative by finite-differencing stencil
 *
 *          The second derivative of the data specified through input to the function is computed using the second order,
 *          central finite-differencing stencil and written into an output array passed to the function.
 *
 * \param   inputMat is the 3D blitz array of which the derivatives have to be computed.
 * \param   outputMat is the 3D blitz array into which the calculated values of derivative are written
 * \param   core is the blitz RectDomain object which offers a view of the core of the matrix
 ********************************************************************************************************************************************
 */
void differ::D2D(blitz::Array<double, 3> inputMat, blitz::Array<double, 3> outputMat, blitz::RectDomain<3> core) {
    outputMat(core) = ((inputMat(shift(core, 1)) - 2.0*inputMat(core) + inputMat(shift(core, -1)))/(h*h));
}
