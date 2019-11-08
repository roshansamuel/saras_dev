#ifndef INITIAL_H
#define INITIAL_H

#include <blitz/array.h>

#include "vfield.h"
#include "grid.h"

class initial {
    private:
        const grid &mesh;

    public:
        initial(const grid &mesh);

        void initTGV(vfield &uField);
        void initRandom(vfield &uField);
        void initSinusoidal(vfield &uField);
};

/**
 ********************************************************************************************************************************************
 *  \class initial initial.h "lib/initial/initial.h"
 *  \brief Contains all the global variables related to the imposing of initial conditions, and functions to impose them
 *
 ********************************************************************************************************************************************
 */

#endif
