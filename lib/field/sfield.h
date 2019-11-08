#ifndef SFIELD_H
#define SFIELD_H

#include "field.h"
#include "derivative.h"

class vfield;       // FORWARD DECLARATION

class sfield {
    private:
        const grid &gridData;
        blitz::Array<double, 3> tempMat;
        
    public:
        field F;
        derivative derS;

        std::string fieldName;

        // The following public arrays for getting interpolated values of variables are available *only if allocDerivatives flag is set to true*
        // Attempting to use these arrays of an sfield with allocDerivatives set to false may give seg-fault!
        blitz::Array<double, 3> interVx, interVy, interVz;
        

        sfield(const grid &gridData, std::string fieldName, const bool allocDerivatives);

        void computeDiff(sfield &H);
        void computeNLin(const vfield &V, sfield &H);

        void gradient(vfield &gradF);

        void syncData();

        sfield& operator += (sfield &a);
        sfield& operator -= (sfield &a);
        sfield& operator *= (double a);

        void operator = (sfield &a);
        void operator = (double a);

        ~sfield() { };
};

/**
 ********************************************************************************************************************************************
 *  \class sfield sfield.h "lib/sfield.h"
 *  \brief Scalar field class to store and operate on scalar fields
 *
 *  The class stores scalar fields in the form of an instance of the field class defined in <B>field.h</B>.
 *  While the class <B>field</B> merely stores data in the form of a blitz array and offers functions to compute derivatives
 *  over a uniform grid, the <B>sfield</B> class adds another layer of functionality along with the <B>grid</B> (<B>grid.h</B>)
 *  class to apply necessary grid transformation metrics and compute derivatives over a non-uniform grid.
 *  The scalar field is also equipped with a gradient operator, which returns a vector field (vfield).
 *  However, this operation is presently restricted to cell-centered scalar fields, i.e., those which are staggered in all
 *  the directions.
 *  Moreover, the \f$ (\mathbf{u}.\nabla)f \f$ operator is also provided as the function <B>computeNLin</B>
 ********************************************************************************************************************************************
 */

#endif
