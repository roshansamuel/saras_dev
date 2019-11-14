#ifndef FORCE_H
#define FORCE_H

#include <blitz/array.h>

#include "parallel.h"
#include "plainsf.h"
#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"

class sfield;

class force{
    public:
        force(vfield &U, const parser &solParams, parallel &mpiParam);

        void add_SForce(plainsf &Ht);
        void add_VForce(plainvf &Hv);
        void add_VForce(plainvf &Hv, sfield &T);

    private:
        vfield &V;

        const parallel &mpiData;

        const parser &inputParams;

        double Fb, Fr;

        void add_Coriolis(plainvf &Hv);
        void add_RandomForce(plainvf &Hv);
        void add_Buoyancy(plainvf &Hv, sfield &T);
};

/**
 ********************************************************************************************************************************************
 *  \class force force.h "lib/force/force.h"
 *  \brief Contains all the global variables related to the imposing of forcing, and associated functions
 *
 ********************************************************************************************************************************************
 */

#endif
