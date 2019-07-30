#ifndef DIFFER_H
#define DIFFER_H

#include <blitz/array.h>

class differ {
    private:
        /** Integer value used to define the axis along which differentation will be done: 0 -> X, 1 -> Y, 2 -> Z */
        const int dim;

        /** Double precision value of grid spacing used in the denominator of finite difference formulae */
        const double h;

    public:
        differ(int dim, double h);

        blitz::RectDomain<3> shift(blitz::RectDomain<3> core, int steps);

        void D1D(blitz::Array<double, 3> inputMat, blitz::Array<double, 3> outputMat, blitz::RectDomain<3> core);
        void D2D(blitz::Array<double, 3> inputMat, blitz::Array<double, 3> outputMat, blitz::RectDomain<3> core);
};

/**
 ********************************************************************************************************************************************
 *  \class differ differ.h "lib/differ.h"
 *  \brief Differ class implements the finite-differencing operations for computing derivatives
 *
 *  The class uses just the grid-spacing along any given direction and an integer value to determine the direction along
 *  which an instance of the class must compute derivatives.
 *  Consequently, the derivatived computed from within this class are for a uniform grid, i.e., it operates on the transformed
 *  coordinates of the computational domain.
 *  The output from this class is used along with a blitz array to form the \ref field data-structure.
 *  To operate on a non-uniform grid, the grid transformation derivatives are used along with the \ref field structure to
 *  form the \ref sfield structure and subquently \ref vfield.
 ********************************************************************************************************************************************
 */

#endif
