#ifndef PLAINSF_H
#define PLAINSF_H

#include "sfield.h"
#include "grid.h"

class plainsf {
    private:
        const grid &gridData;

    public:
        blitz::Array<double, 3> F;

        plainsf(const grid &gridData, const sfield &refF);

        inline double fieldMax();

        plainsf& operator += (plainsf &a);
        plainsf& operator -= (plainsf &a);

        plainsf& operator += (sfield &a);
        plainsf& operator -= (sfield &a);

        plainsf& operator *= (double a);

        void operator = (plainsf &a);
        void operator = (sfield &a);

        void operator = (double a);

        ~plainsf() { };
};

/**
 ********************************************************************************************************************************************
 *  \class plainsf plainsf.h "lib/plainsf.h"
 *  \brief Plain scalar field class to store simple scalar fields with no differentiation or interpolation
 *
 *  The class stores scalar fields in the form of a Blitz array
 ********************************************************************************************************************************************
 */

#endif
