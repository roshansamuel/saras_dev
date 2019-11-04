#ifndef FORCE_H
#define FORCE_H

#include <blitz/array.h>
#include <string>

#include "parallel.h"
#include "plainsf.h"
#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"

class sfield;

class force{
    private:
        double Fb, Fr;
        void add_Buoyancy(plainvf &Hv, sfield &T);
        void add_Coriolis(plainvf &Hv);
        void add_RandomForce(plainvf &Hv);

    public:
        vfield V;
        //sfield T;
        const parser &inputParams;
        parallel &mpiData;

        force(vfield &U, const parser &solParams, parallel &mpiParam);
        //force(vfield &U, sfield &S, const parser &solParams, parallel &mpiParam);

        void add_VForce(plainvf &Hv);
        void add_VForce(plainvf &Hv, sfield &T);
        void add_SForce(plainsf &Ht);

        ~force();
};



#endif
        

