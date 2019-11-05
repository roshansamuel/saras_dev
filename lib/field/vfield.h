#ifndef VFIELD_H
#define VFIELD_H

#include "field.h"

// Forward declarations of relevant classes
class plainsf;
class plainvf;
class sfield;

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

        vfield(const grid &gridData, std::string fieldName);

        void computeDiff(plainvf &H);
        void computeNLin(const vfield &V, plainvf &H);
        void computeTStp(double &dt_out);

        void divergence(plainsf &divV, const sfield &P);

        void syncData();

        vfield& operator += (plainvf &a);
        vfield& operator -= (plainvf &a);

        vfield& operator += (vfield &a);
        vfield& operator -= (vfield &a);

        vfield& operator *= (double a);

        void operator = (plainvf &a);
        void operator = (vfield &a);
        void operator = (double a);

        ~vfield() { };
};

/**
 ********************************************************************************************************************************************
 *  \class vfield vfield.h "lib/vfield.h"
 *  \brief Vector field class to store and operate on vector fields
 *
 *  The class stores vector fields in the form of three instances of the field class defined in <B>field.h</B>.
 *  The vector field is stored in such a way that the components are face-centered, with:
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
