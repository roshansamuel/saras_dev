#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <blitz/array.h>

#include "sfield.h"
#include "grid.h"

class boundary {
    private:
        const grid &mesh;

        sfield &dField;

        bool nonHgBC;

        blitz::Array<double, 1> x, y, z;
        blitz::Array<double, 1> xGlo, yGlo, zGlo;

        void setXYZ();

    public:
        blitz::Array<bool, 3> wallMask;
        blitz::Array<double, 3> wallData;

        /*****************************************************************************************************************************************************/

        boundary(const grid &mesh, sfield &inField);

        void createPatch(int wallNum);

        void imposeBC();
};

/**
 ********************************************************************************************************************************************
 *  \class boundary boundary.h "lib/boundary/boundary.h"
 *  \brief Contains all the global variables related to the imposing of boundary conditions, and functions to impose BCs
 *
 ********************************************************************************************************************************************
 */

#endif
