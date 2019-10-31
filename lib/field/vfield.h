#ifndef VFIELD_H
#define VFIELD_H

#include "field.h"
#include "parallel.h"
#include "grid.h"
#include <math.h>

class sfield;       // FORWARD DECLARATION

class vfield {
    private:
        const grid &gridData;

    public:
        field Vx, Vy, Vz;

        std::string fieldName;

        // The following public arrays for getting interpolated values of variables are available *only if allocDerivatives flag is set to true*
        // Attempting to use these arrays of a vfield with allocDerivatives set to false may give seg-fault!
        blitz::Array<double, 3> interVx2Vx, interVy2Vx, interVz2Vx;
        blitz::Array<double, 3> interVx2Vy, interVy2Vy, interVz2Vy;
        blitz::Array<double, 3> interVx2Vz, interVy2Vz, interVz2Vz;
        blitz::Array<double, 3> interPc2Vz;

        vfield(const grid &gridData, std::string fieldName, const bool allocDerivatives);

        void computeDiff(vfield &H);
        void computeTStp(double &dt_out, double Courant_no);
        void computeNLin(const vfield &V, vfield &H);

        void divergence(sfield &divV);

        void syncData();

        vfield& operator += (vfield &a);
        vfield& operator -= (vfield &a);
        vfield& operator *= (double a);

        void operator = (vfield &a);
        void operator = (double a);

        ~vfield() { };
};

/**
 ********************************************************************************************************************************************
 *  \class vfield vfield.h "lib/vfield.h"
 *  \brief Vector field class to store and operate on vector fields
 *
 *  The class stores vector fields in the form of three instances of the sfield class defined in <B>sfield.h</B>.
 *  The vector field is stored in such a way that the components are face-centered scalar fields, with:
 *      - x-component located at the face centers along the yz-plane
 *      - y-component located at the face centers along the zx-plane
 *      - z-component located at the face centers along the xy-plane
 *
 *  The vector field is also equipped with a divergence operator, which returns a scalar field (sfield).
 *  However, this operation returns only cell-centered scalar field as output as most scalar fields are stored at
 *  cell centers.
 *  Moreover, the \f$ (\mathbf{u}.\nabla)\mathbf{v} \f$ operator is also provided as the function <B>computeNLin</B>
 ********************************************************************************************************************************************
 */

#endif
