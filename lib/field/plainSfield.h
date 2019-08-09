#ifndef PSFIELD_H
#define PSFIELD_H

#include "field.h"

class pvfield;       // FORWARD DECLARATION

class psfield {
    
    private:
        const grid &gridData;
        
    public:
        field F;

        std::string fieldName;

        
        psfield(const grid &gridData, std::string fieldName);

        void syncData();

        psfield& operator += (psfield &a);
        psfield& operator -= (psfield &a);
        psfield& operator *= (double a);

        void operator = (psfield &a);
        void operator = (double a);

        ~psfield() { };
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
