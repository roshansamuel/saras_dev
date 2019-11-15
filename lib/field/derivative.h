#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <blitz/array.h>
#include <blitz/array/stencil-et.h>
#include <blitz/array/stencilops.h>
#include <string>

#include "field.h"
#include "grid.h"

class derivative {
    private: 
        const grid &gridData;

        const field &F;

        double invDelx, invDely, invDelz;

        blitz::firstIndex i;
        blitz::secondIndex j;
        blitz::thirdIndex k;    

        blitz::Range fullRange;

        blitz::Array<double, 1> x_Metric, y_Metric, z_Metric;
        blitz::Array<double, 1> xxMetric, yyMetric, zzMetric;
        blitz::Array<double, 1> x2Metric, y2Metric, z2Metric;

        blitz::Array<double, 3> tempMat;

    public:
        derivative(const grid &gridData, const field &F);

        void calcDerivative1_x(blitz::Array<double, 3> outputMat);
        void calcDerivative1_y(blitz::Array<double, 3> outputMat);
        void calcDerivative1_z(blitz::Array<double, 3> outputMat);

        void calcDerivative2xx(blitz::Array<double, 3> outputMat);
        void calcDerivative2yy(blitz::Array<double, 3> outputMat);
        void calcDerivative2zz(blitz::Array<double, 3> outputMat);
};

/**
 ********************************************************************************************************************************************
 *  \class derivative derivative.h "lib/derivative.h"
 *  \brief Derivative class to perform finite difference operations on the data stored in field
 *
 *  It contains functions to perform the finite difference operations with constant grid spacing.
 *  For many classes of this solver, empty destructors are removed. Refer reference [5] of README for more details.
 ********************************************************************************************************************************************
 */

#endif
