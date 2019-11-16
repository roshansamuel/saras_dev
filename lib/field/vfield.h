#ifndef VFIELD_H
#define VFIELD_H

#include "field.h"
#include "derivative.h"

// Forward declarations of relevant classes
class plainsf;
class plainvf;
class sfield;

class vfield {
    private:
        const grid &gridData;

        blitz::Array<double, 3> derivTempX, derivTempY, derivTempZ;

    public:
        field Vx, Vy, Vz;

        derivative derVx, derVy, derVz;

        std::string fieldName;

        blitz::Array<double, 3> interTempX, interTempY, interTempZ;

        vfield(const grid &gridData, std::string fieldName);

        void computeDiff(plainvf &H);
        void computeTStp(double &dt_out);
        void computeNLin(const vfield &V, plainvf &H);

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
