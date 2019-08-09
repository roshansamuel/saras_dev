#ifndef PVFIELD_H
#define PVFIELD_H

#include "field.h"
#include "parser.h"
#include "parallel.h"
#include "grid.h"
#include <math.h>

class sfield;       // FORWARD DECLARATION

class pvfield {
    private:
        const grid &gridData;

    public:
        field Vx, Vy, Vz;

        std::string fieldName;

        
        pvfield(const grid &gridData, std::string fieldName);

        void syncData();

        pvfield& operator += (pvfield &a);
        pvfield& operator -= (pvfield &a);
        pvfield& operator *= (double a);

        void operator = (pvfield &a);
        void operator = (double a);

        ~pvfield() { };
};

/**
 ********************************************************************************************************************************************
 *  \class plainVfield plainVfield.h "lib/plainVfield.h"
 *  \brief Vector field class to store and operate vector fields
 *
 *  The class stores vector fields in the form of three instances of the plainSfield class defined in <B>plainSfield.h</B>.
 *  The vector field is stored in such a way that the components are face-centered scalar fields, with:
 *      - x-component located at the face centers along the yz-plane
 *      - y-component located at the face centers along the zx-plane
 *      - z-component located at the face centers along the xy-plane
 *
 *  
 ********************************************************************************************************************************************
 */

#endif
