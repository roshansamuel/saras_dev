#ifndef FORCE_H
#define FORCE_H

#include <blitz/array.h>
#include <string>

#include "parallel.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"

class sfield;

class force{
    private:
        double Fb, Fr;
        void add_Buoyancy(vfield &Hv, sfield &T);
        void add_Coriolis(vfield &Hv);
        void add_RandomForce(vfield &Hv);

    public:
        //blitz::Array<double, 3> Force_x, Force_y, Force_z;
        vfield V;
        //sfield T;
        const parser &inputParams;
        parallel &mpiData;

        force(vfield &U, const parser &solParams, parallel &mpiParam);
        //force(vfield &U, sfield &S, const parser &solParams, parallel &mpiParam);

        void add_VForce(vfield &Hv);
        void add_VForce(vfield &Hv, sfield &T);
        void add_SForce(sfield &Ht);

        ~force();
};



#endif
        

